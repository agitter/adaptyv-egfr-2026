"""Matched-component MM/GBSA checks of interface histidine protonation.

Amber99SB/OBC, fixed heavy-atom coordinates. Energies are complex minus both
partners WITH IDENTICAL protonation and coordinates. Different-proton-count
absolute force-field energies are never compared as pH-dependent affinities.
"""
from pathlib import Path
import sys,os,json,time,copy,random
import numpy as np
from mm_refine import *
from classical_design import LIB,conformer,pair_score,atom_data,O,charge_center
from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1


def complete_target(rows):
    rows=copy.deepcopy(rows)
    for i,r in enumerate(rows):
        if set(LIB[r['aa']]['names'])<=set(r['atoms']):continue
        bb=np.array([r['atoms'][n] for n in BB]+[r['atoms'].get('CB',r['atoms']['CA'])])
        env=np.concatenate([atom_data(a['aa'],a['atoms']) for j,a in enumerate(rows) if j!=i and np.linalg.norm(np.array(a['atoms']['CA'])-bb[1])<15])
        options=[]
        for k in range(len(LIB[r['aa']]['xyz'])):
            a=conformer(r['aa'],k,bb);x=atom_data(r['aa'],a,LIB[r['aa']]['names']);rep,vw,hb=pair_score(x,env)
            fit=sum(np.sum((a[n]-p)**2) for n,p in r['atoms'].items() if n not in BB and n in a)
            options.append((rep+fit+vw+hb,a))
        _,atoms=min(options,key=lambda a:a[0]);r['atoms']={k:np.array(v).tolist() for k,v in atoms.items()}
    return rows


def energy(top,pos):
    system=FF.createSystem(top,nonbondedMethod=app.NoCutoff,constraints=None,removeCMMotion=False)
    integ=mm.VerletIntegrator(.001*u.picoseconds);ctx=mm.Context(system,integ,mm.Platform.getPlatformByName('CPU'),{'Threads':'1'});ctx.setPositions(pos)
    val=ctx.getState(getEnergy=True).getPotentialEnergy().value_in_unit(u.kilocalorie_per_mole)
    del ctx,integ,system
    return float(val)


def mmgbsa(heavy,variants):
    random.seed(19982003);np.random.seed(19982003)
    mod=app.Modeller(heavy.topology,heavy.positions)
    used=mod.addHydrogens(FF,pH=7.4,variants=variants,platform=mm.Platform.getPlatformByName('CPU'))
    total=energy(mod.topology,mod.positions);components={}
    for chain in ['A','B']:
        part=app.Modeller(mod.topology,mod.positions);part.delete([a for a in part.topology.atoms() if a.residue.chain.id!=chain]);components[chain]=energy(part.topology,part.positions)
    return {'delta_kcal_proxy':total-components['A']-components['B'],'complex_kcal':total,'target_kcal':components['A'],'binder_kcal':components['B'],'variants':used}


def main(candidate='P00024_07',species='human6ARU',single_only=True):
    out=S/'intermediate/mm_binding';out.mkdir(exist_ok=True)
    d=json.loads((S/'intermediate/designs'/(candidate+'.json')).read_text());rows=d['structure']
    relaxed=S/'intermediate/refined'/(candidate+'_relaxed.pdb')
    if relaxed.exists():
        pdb=PDBParser(QUIET=True).get_structure('binder',str(relaxed));atoms=[{a.name:a.coord.tolist() for a in r if a.element!='H'} for r in next(pdb[0].get_chains()) if 'CA' in r]
        assert len(atoms)==len(rows)
        rows=copy.deepcopy(rows)
        for r,a in zip(rows,atoms):r['atoms']=a
    target=json.loads((O/(species+'_aligned.json')).read_text());target=[r for r in target if 334<=r['human_pos']<=505]
    assert len(target)==172,'Need explicit repair of target chain gaps'
    target=complete_target(target)
    bp=out/(candidate+'_binder.pdb');tp=out/(candidate+'_'+species+'_target.pdb');write_pdb(rows,bp,'B');write_pdb(target,tp,'A')
    t=app.PDBFile(str(tp));b=app.PDBFile(str(bp));heavy=app.Modeller(t.topology,t.positions);heavy.add(b.topology,b.positions);heavy.topology.createDisulfideBonds(heavy.positions)
    joint=out/(candidate+'_'+species+'_complex_relaxed.pdb')
    if joint.exists():
        jp=app.PDBFile(str(joint));heavy=app.Modeller(jp.topology,jp.positions)
        heavy.delete([a for a in heavy.topology.atoms() if a.element==app.element.hydrogen])
    variants=[None]*heavy.topology.getNumResidues();hist=[]
    xb=np.array([v for r in rows for v in r['atoms'].values()]);xt=np.array([v for r in target for v in r['atoms'].values()])
    for i,r in enumerate(target+rows):
        if r['aa']!='H':continue
        variants[i]='HIE';center=charge_center('H',r['atoms']);other=xb if i<len(target) else xt
        distance=np.linalg.norm(other-center,axis=1).min()
        if distance<8.0:hist.append({'topology_index':i,'partner':'target' if i<len(target) else 'binder','residue_position':r.get('human_pos',r.get('position')),'nearest_partner_distance':float(distance)})
    tic=time.time();result={'candidate_id':candidate,'species':species,'binding_geometry':'jointly relaxed complex' if joint.exists() else ('relaxed binder only' if relaxed.exists() else 'unrelaxed designed binder'),'target_crop':[334,505],'histidines':hist,'single_states':[],'interpretation':'MM/GBSA diagnostic, not measured affinity or complete proton-linked binding free energy'}
    base=mmgbsa(heavy,variants);result['neutral_reference']=base
    for h in hist:
        index=h['topology_index'];v=variants.copy();v[index]='HIP';e=mmgbsa(heavy,v)
        result['single_states'].append({**h,'protonated':e,'delta_delta_binding_kcal_proxy':e['delta_kcal_proxy']-base['delta_kcal_proxy']})
        print(candidate,species,h['partner'],h['residue_position'],'protonation binding shift',round(e['delta_kcal_proxy']-base['delta_kcal_proxy'],3),flush=True)
    result['seconds']=time.time()-tic
    (out/(candidate+'_'+species+('_joint_probe.json' if joint.exists() else '_probe.json'))).write_text(json.dumps(result,indent=2));print('MMGBSA completed',candidate,species,round(time.time()-tic,2),flush=True)
if __name__=='__main__':main(*sys.argv[1:])
