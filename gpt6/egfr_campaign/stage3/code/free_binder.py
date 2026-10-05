"""Local unbound relaxation / optional deterministic Verlet stability check.

This assesses reorganization and obvious local instability, NOT de novo folding
or an equilibrium binding free energy. All algorithms predate 2011. Stochastic
initial velocities, when requested, use explicit MT19937, not an engine RNG.
"""
from pathlib import Path
import sys,json,time,math
import numpy as np
R=Path('/mnt/data/egfr_campaign');S3=R/'stage3';sys.path[:0]=[str(R/'stage2/code'),str(R/'stage2/runtime/python'),str(R/'code')]
from mm_refine import *


def align_metrics(top,pos,initial,mask):
    xyz=np.asarray(pos.value_in_unit(u.angstrom));x0=np.asarray(initial.value_in_unit(u.angstrom));ids=np.array([a.index for a in top.atoms() if a.name=='CA']);x=x0[ids];y=xyz[ids];fm=~np.asarray(mask,dtype=bool)
    center0=x[fm].mean(0);center=y[fm].mean(0);U,ss,V=np.linalg.svd((y[fm]-center).T@(x[fm]-center0));rot=U@V
    if np.linalg.det(rot)<0:U[:,-1]*=-1;rot=U@V
    fit=(y-center)@rot+center0;ds=np.linalg.norm(fit-x,axis=1);starts=np.flatnonzero(np.diff(np.r_[False,mask,False]).astype(int));cdr=[]
    for a,b in zip(starts[::2],starts[1::2]):cdr.append(float(np.sqrt(np.mean(ds[a:b]**2))))
    return {'framework_CA_RMSD_A':float(np.sqrt(np.mean(ds[fm]**2))),'loop_CA_RMSD_A':float(np.sqrt(np.mean(ds[~fm]**2))) if (~fm).sum() else 0.,'CDR_CA_RMSD_A':cdr,'maximum_CA_displacement_A':float(ds.max())}


def main(cid,steps=0):
    steps=int(steps);out=S3/'intermediate/unbound';out.mkdir(exist_ok=True);d=json.loads((S3/'intermediate/designs'/(cid+'.json')).read_text());source=S3/'intermediate/refined'/(cid+'_human6ARU_relaxed.pdb');pdb=app.PDBFile(str(source));mod=app.Modeller(pdb.topology,pdb.positions);mod.delete([r for r in mod.topology.residues() if r.chain.id!='B']);mod.topology.createDisulfideBonds(mod.positions)
    initial=mod.positions;ref=chirality_groups(mod.topology,initial);mask=[r['loop'] for r in d['structure']];assert len(mask)==mod.topology.getNumResidues()
    system=FF.createSystem(mod.topology,nonbondedMethod=app.NoCutoff,constraints=app.HBonds,removeCMMotion=False)
    for f in system.getForces():f.setForceGroup(0)
    restraint=mm.CustomTorsionForce('k*(1-cos(theta-theta0))');restraint.addPerTorsionParameter('theta0');restraint.addGlobalParameter('k',1000.);restraint.setForceGroup(1)
    for kind,res,ids,theta in ref:restraint.addTorsion(*ids,[theta])
    system.addForce(restraint);integ=mm.VerletIntegrator(.001*u.picoseconds);ctx=mm.Context(system,integ,mm.Platform.getPlatformByName('CPU'),{'Threads':'1'});ctx.setPositions(initial);start=time.time();e0=float(ctx.getState(getEnergy=True,groups=1).getPotentialEnergy().value_in_unit(u.kilocalorie_per_mole));mm.LocalEnergyMinimizer.minimize(ctx,10.,1000);state=ctx.getState(getPositions=True,getEnergy=True);end=state.getPositions();e1=float(ctx.getState(getEnergy=True,groups=1).getPotentialEnergy().value_in_unit(u.kilocalorie_per_mole));audit=geometry_audit(mod.topology,end,ref);metrics=align_metrics(mod.topology,end,initial,mask)
    with open(out/(cid+'_free_relaxed.pdb'),'w') as f:app.PDBFile.writeFile(mod.topology,end,f)
    result={'candidate_id':cid,'source':str(source),'minimization_iterations_limit':1000,'position_restraints':False,'chirality_restraints':True,'energy_excludes_chirality_restraints':True,'before_physical_potential_kcal':e0,'after_physical_potential_kcal':e1,'potential_relaxation_kcal':e0-e1,'geometry_audit':audit,'displacements':metrics,'interpretation':'Local unbound reorganization diagnostic; not folding validation or binding entropy.'}
    if steps:
        rng=np.random.Generator(np.random.MT19937(198019672003));mass=np.array([system.getParticleMass(i).value_in_unit(u.dalton) for i in range(system.getNumParticles())]);vel=rng.normal(size=(len(mass),3))*np.sqrt(.008314462618*300./mass)[:,None];vel-=(vel*mass[:,None]).sum(0)/mass.sum();ctx.setVelocities(vel*u.nanometer/u.picosecond);ctx.applyVelocityConstraints(1e-6);frames=[];all_e=[]
        for step in range(0,steps+1,1000):
            if step:integ.step(min(1000,steps-step+1000))
            st=ctx.getState(getPositions=True,getEnergy=True);potential=float(st.getPotentialEnergy().value_in_unit(u.kilojoule_per_mole));kin=float(st.getKineticEnergy().value_in_unit(u.kilojoule_per_mole));all_e.append(potential+kin);a=geometry_audit(mod.topology,st.getPositions(),ref);frames.append({'step':step,'potential_kj':potential,'kinetic_kj':kin,'total_kj':potential+kin,'geometry_pass':a['pass'],**align_metrics(mod.topology,st.getPositions(),initial,mask)})
        with open(out/(cid+'_free_NVE_final.pdb'),'w') as f:app.PDBFile.writeFile(mod.topology,st.getPositions(),f)
        result['NVE']={'steps':steps,'timestep_ps':.001,'initial_temperature_K':300.,'rng':'MT19937','thermostat':None,'algorithm':'Verlet 1967 with HBond constraints; no modern stochastic integrator','frames':frames,'total_energy_range_kj':float(max(all_e)-min(all_e)),'qualification':'Very short local trajectory, not convergence or a fold prediction.'}
    result['seconds']=time.time()-start;(out/(cid+'_free_audit.json')).write_text(json.dumps(result,indent=2));print(cid,'unbound',audit['pass'],metrics,'potential relaxation',e0-e1,'seconds',round(result['seconds'],1),flush=True);del ctx,integ
if __name__=='__main__':main(*sys.argv[1:])
