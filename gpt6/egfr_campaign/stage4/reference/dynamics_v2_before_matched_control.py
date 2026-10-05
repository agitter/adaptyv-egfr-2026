"""Matched classical thermal stress test, not folding or equilibrium validation.
Amber99SB(2006)/OBC(2004), Verlet(1967), constrained dynamics, Andersen(1980)
style collisions of H-bond constraint groups. Random numbers: MT19937(1998).
"""
from pathlib import Path
import sys,os,json,time,hashlib,math,random
import numpy as np
R=Path('/mnt/data/egfr_campaign');S=R/'stage4'
sys.path[:0]=[str(R/'stage3/code'),str(R/'stage2/code'),str(R/'code'),str(R/'stage2/runtime/python')]
from mm_refine import FF,app,mm,u,chirality_groups,geometry_audit,write_pdb
from free_binder import align_metrics
from classical_design import LIB,BB,conformer,atom_data,pair_score
from Bio.SeqUtils import seq1

def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def prepare(cid):
    random.seed(19801998);np.random.seed(19801998)
    if cid=='NATIVE_3EAK':
        source=R/'inputs/transfer1/legacy_fetch/3EAK.pdb';p=app.PDBFile(str(source));x=np.asarray(p.positions.value_in_unit(u.angstrom));rows=[]
        for res in p.topology.residues():
            if res.chain.id!='A' or res.name in ['HOH','SO4']:continue
            aa=seq1(res.name);rows.append({'aa':aa,'atoms':{a.name:x[a.index].tolist() for a in res.atoms() if a.element!=app.element.hydrogen}})
        assert len(rows)==127
        repairs=[]
        for i,row in enumerate(rows):
            aa=row['aa'];missing=set(LIB[aa]['names'])-set(row['atoms'])
            if not missing:continue
            env=np.concatenate([atom_data(r['aa'],r['atoms']) for j,r in enumerate(rows) if abs(i-j)>1]);bb=np.array([row['atoms'][n] for n in BB]);opts=[]
            for k in range(len(LIB[aa]['xyz'])):
                a=conformer(aa,k,bb);v=pair_score(atom_data(aa,a,LIB[aa]['names']),env)
                observed=sum(np.sum((np.array(row['atoms'][n])-a[n])**2) for n in row['atoms'] if n in a and n not in BB)
                opts.append((sum(v)+observed,k,a))
            val,k,a=min(opts,key=lambda z:z[0]);row['atoms']={n:np.asarray(v).tolist() for n,v in a.items()};repairs.append({'position':i+1,'rotamer':k,'missing':sorted(missing)})
        dest=S/'intermediate/dynamics/native_heavy.pdb';write_pdb(rows,dest,'B');p=app.PDBFile(str(dest));mod=app.Modeller(p.topology,p.positions)
        intervals=[[25,38],[52,67],[99,117]];mask=[any(a<=i<b for a,b in intervals) for i in range(127)]
        meta={'role':'native experimental fold control ONLY; original CDRs not used as design seeds','repairs':repairs,'cdr_mask_qualification':'native source residue regions26-38,53-67,100-117; loop lengths differ from designs'}
    else:
        root=next(r for r in [S,R/'stage3'] if (r/'intermediate/refined'/f'{cid}_human6ARU_relaxed.pdb').exists());source=root/'intermediate/refined'/f'{cid}_human6ARU_relaxed.pdb';d=json.loads((root/'intermediate/designs'/f'{cid}.json').read_text());p=app.PDBFile(str(source));mod=app.Modeller(p.topology,p.positions);mod.delete([r for r in mod.topology.residues() if r.chain.id!='B']);mask=[r['loop'] for r in d['structure']];meta={'role':'untested de novo computational design','sequence':d['sequence']}
    mod.topology.createDisulfideBonds(mod.positions);mod.delete([a for a in mod.topology.atoms() if a.element==app.element.hydrogen]);mod.addHydrogens(FF,pH=7.4,platform=mm.Platform.getPlatformByName('CPU'))
    assert len(mask)==mod.topology.getNumResidues()
    return mod,mask,{'source':str(source),'source_sha256':digest(source),**meta}

def constraint_groups(system):
    groups={i:{i} for i in range(system.getNumParticles())};labels=list(range(system.getNumParticles()))
    for k in range(system.getNumConstraints()):
        i,j,d=system.getConstraintParameters(k);a,b=labels[i],labels[j]
        if a!=b:
            for n in groups[b]:labels[n]=a
            groups[a]|=groups.pop(b)
    return [sorted(v) for v in groups.values()]

def velocity_verlet(dt_ps):
    """Synchronized velocity Verlet with RATTLE-style constraint corrections.
    Scientific equations predate2011; CustomIntegrator is only the implementation.
    """
    it=mm.CustomIntegrator(dt_ps*u.picosecond)
    it.addPerDofVariable('x_free',0.)
    it.addUpdateContextState()
    it.addComputePerDof('v','v+0.5*dt*f/m')
    it.addComputePerDof('x','x+dt*v')
    it.addComputePerDof('x_free','x')
    it.addConstrainPositions()
    it.addComputePerDof('v','v+0.5*dt*f/m+(x-x_free)/dt')
    it.addConstrainVelocities()
    it.setConstraintTolerance(1e-6)
    return it

def run(cid,seed=198001,production_ps=4.,warm_ps=2.,threads=2,chirality=True):
    start=time.time();out=S/'intermediate/dynamics';out.mkdir(parents=True,exist_ok=True);tag=f'{cid}_v2_s{seed}';dest=out/(tag+'.json')
    if dest.exists():raise FileExistsError('Refuse overwrite '+str(dest))
    rng=np.random.Generator(np.random.MT19937(int(seed)));mod,mask,meta=prepare(cid);initial=mod.positions;ref=chirality_groups(mod.topology,initial)
    system=FF.createSystem(mod.topology,nonbondedMethod=app.NoCutoff,constraints=app.HBonds,removeCMMotion=False)
    if chirality:
        force=mm.CustomTorsionForce('kc*(1-cos(theta-theta0))');force.addGlobalParameter('kc',1000.);force.addPerTorsionParameter('theta0')
        for kind,res,ids,theta in ref:force.addTorsion(*ids,[theta])
        force.setForceGroup(1);system.addForce(force)
    restraint=mm.CustomExternalForce('ramp*kpos*((x-x0)^2+(y-y0)^2+(z-z0)^2)');restraint.addGlobalParameter('ramp',1.)
    for par in ['kpos','x0','y0','z0']:restraint.addPerParticleParameter(par)
    pos=np.asarray(initial.value_in_unit(u.nanometer))
    for a in mod.topology.atoms():
        if a.name in BB:restraint.addParticle(a.index,[500.,*pos[a.index]])
    restraint.setForceGroup(2);system.addForce(restraint)
    dt=.002;integrator=velocity_verlet(dt);integrator.setConstraintTolerance(1e-6);ctx=mm.Context(system,integrator,mm.Platform.getPlatformByName('CPU'),{'Threads':str(threads)});ctx.setPositions(initial);mm.LocalEnergyMinimizer.minimize(ctx,10.,500)
    groups=constraint_groups(system);dof=3*system.getNumParticles()-system.getNumConstraints()-3;mass=np.array([system.getParticleMass(i).value_in_unit(u.dalton) for i in range(system.getNumParticles())]);sigma=np.sqrt(.008314462618*300/mass)
    vel=rng.normal(size=(len(mass),3))*sigma[:,None];vel-=(vel*mass[:,None]).sum(0)/mass.sum();ctx.setVelocities(vel*u.nanometer/u.picosecond);ctx.applyVelocityConstraints(1e-6)
    interval_steps=25;prob=1-math.exp(-5.*dt*interval_steps);n_warm=int(round(float(warm_ps)/dt));n_prod=int(round(float(production_ps)/dt));assert n_warm%interval_steps==0 and n_prod%interval_steps==0
    frames=[];trace=[];ca_ids=[a.index for a in mod.topology.atoms() if a.name=='CA']
    for step in range(0,n_warm+n_prod+1,interval_steps):
        warm=step<n_warm;ramp=max(0.,1-step/max(n_warm,1)) if warm else 0.;ctx.setParameter('ramp',ramp)
        if step%250==0 or step==n_warm+n_prod:
            st=ctx.getState(getPositions=True,getEnergy=True);pp=st.getPositions();ke=float(st.getKineticEnergy().value_in_unit(u.kilojoule_per_mole));metric=align_metrics(mod.topology,pp,initial,mask)
            frame={'step':step,'time_ps':step*dt,'phase':'warmup' if warm else 'production','ramp':ramp,'temperature_K':2*ke/(dof*.008314462618),'potential_kj':float(st.getPotentialEnergy().value_in_unit(u.kilojoule_per_mole)),**metric};frames.append(frame);trace.append(np.asarray(pp.value_in_unit(u.angstrom))[ca_ids]);(out/(tag+'_progress.json')).write_text(json.dumps({'candidate_id':cid,'frames':frames,'elapsed_seconds':time.time()-start},indent=2))
            if not np.isfinite(frame['potential_kj']):raise ValueError('nonfinite energy')
            print(cid,seed,round(step*dt,2),frame['phase'],'T',round(frame['temperature_K']),'CDRs',[round(x,2) for x in metric['CDR_CA_RMSD_A']],flush=True)
        if step==n_warm+n_prod:break
        velocities=np.asarray(ctx.getState(getVelocities=True).getVelocities(asNumpy=True).value_in_unit(u.nanometer/u.picosecond));sel=rng.random(len(groups))<prob;ids=[i for k,g in enumerate(groups) if sel[k] for i in g]
        if ids:velocities[ids]=rng.normal(size=(len(ids),3))*sigma[ids,None]
        velocities-=(velocities*mass[:,None]).sum(0)/mass.sum();ctx.setVelocities(velocities*u.nanometer/u.picosecond);ctx.applyVelocityConstraints(1e-6);integrator.step(interval_steps)
    final=ctx.getState(getPositions=True).getPositions();audit=geometry_audit(mod.topology,final,ref)
    with open(out/(tag+'_final.pdb'),'w') as f:app.PDBFile.writeFile(mod.topology,final,f)
    np.savez_compressed(out/(tag+'_CA_trace.npz'),coords_A=np.array(trace),times_ps=np.array([f['time_ps'] for f in frames]))
    prod=[f for f in frames if f['phase']=='production'];result={'candidate_id':cid,'seed':seed,'rng':'MT19937','warmup_ps':warm_ps,'production_ps':production_ps,'timestep_ps':dt,'protocol':'v2:Andersen-style constraint-group collisions every0.05ps,5/ps;synchronized velocity-Verlet/RATTLE;300K;NoCutoff','forcefield':'Amber99SB2006/OBC2004','chirality_restraints':chirality,'position_restraints_production':False,'metadata':meta,'frames':frames,'summary':{'production_mean_temperature_K':float(np.mean([f['temperature_K'] for f in prod])),'production_mean_loop_RMSD_A':float(np.mean([f['loop_CA_RMSD_A'] for f in prod])),'final_CDR_RMSD_A':prod[-1]['CDR_CA_RMSD_A'],'final_framework_RMSD_A':prod[-1]['framework_CA_RMSD_A'],'production_mean_CDR_RMSD_A':np.mean([f['CDR_CA_RMSD_A'] for f in prod],0).tolist()},'final_geometry_audit':audit,'seconds':time.time()-start,'limitation':'Very short nonconverged thermal stress test, not equilibrium folding stability, physical kinetics or binding free energy.'}
    dest.write_text(json.dumps(result,indent=2));print('DONE',tag,result['summary'],'seconds',result['seconds'],flush=True);del ctx,integrator
if __name__=='__main__':run(sys.argv[1],int(sys.argv[2]) if len(sys.argv)>2 else 198001,float(sys.argv[3]) if len(sys.argv)>3 else 4.,float(sys.argv[4]) if len(sys.argv)>4 else 2.)
