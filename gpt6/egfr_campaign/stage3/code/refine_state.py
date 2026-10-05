"""Refine new candidate interfaces, audit them, then enumerate proton states."""
from pathlib import Path
import sys,json,time,copy
R=Path('/mnt/data/egfr_campaign');S3=R/'stage3'
sys.path[:0]=[str(S3/'code'),str(R/'stage2/code'),str(R/'stage2/runtime/python'),str(R/'code')]
from explicit_refine import *
from mm_binding_probe import complete_target
from proton_ensemble import probe
from scipy.spatial import cKDTree


def run(candidate,species='human6ARU',iterations=600):
    out=S3/'intermediate/refined';out.mkdir(parents=True,exist_ok=True);name=candidate+'_'+species+'_acid';prefix=out/name
    if Path(str(prefix)+'_proton.json').exists():return {'candidate_id':candidate,'status':'cached'}
    d=json.loads((S3/'intermediate/designs'/(candidate+'.json')).read_text());rows=d['structure'];alltarget=json.loads((R/'intermediate'/(species+'_aligned.json')).read_text());target=complete_target([r for r in alltarget if 334<=r['human_pos']<=505]);assert len(target)==172
    tp=Path(str(prefix)+'_target.pdb');bp=Path(str(prefix)+'_binder.pdb');write_pdb(target,tp,'A');write_pdb(rows,bp,'B')
    a=app.PDBFile(str(tp));b=app.PDBFile(str(bp));mod=app.Modeller(a.topology,a.positions);mod.add(b.topology,b.positions)
    src=Path(str(prefix)+'_initial.pdb');dst=Path(str(prefix)+'_relaxed.pdb')
    with open(src,'w') as f:app.PDBFile.writeFile(mod.topology,mod.positions,f)
    # Human glycan positions also conservatively protect the homologous mouse surface.
    gly=json.loads((R/'stage2/intermediate/target_glycans.json').read_text());mask=[False]*len(target)+[r['loop'] for r in rows]
    audit_path=Path(str(prefix)+'_audit.json')
    if not audit_path.exists():
        result=minimize(src,dst,mask,iterations,[g['xyz'] for g in gly]);result.update({'structural_hypothesis':'acid all-interface-HIP; pH dependence calculated separately', 'candidate_id':candidate,'species':species,'target_crop':[334,505],'joint_refinement':True,'glycan_obstacle_count':len(gly)})
        p=app.PDBFile(str(dst));xyz=np.asarray(p.positions.value_in_unit(u.angstrom));bind=[a.index for a in p.topology.atoms() if a.residue.chain.id=='B' and a.element!=app.element.hydrogen]
        outside=np.array([v for r in alltarget if not 334<=r['human_pos']<=505 for v in r['atoms'].values()]);gxyz=np.array([g['xyz'] for g in gly])
        result['postrefinement_full_receptor_outside_crop_min_A']=float(cKDTree(outside).query(xyz[bind])[0].min());result['postrefinement_glycan_min_A']=float(cKDTree(gxyz).query(xyz[bind])[0].min())
        result['postrefinement_exclusion_pass']=result['postrefinement_full_receptor_outside_crop_min_A']>=2. and result['postrefinement_glycan_min_A']>=2.
        audit_path.write_text(json.dumps(result,indent=2))
    else:result=json.loads(audit_path.read_text())
    print(name,'refinement',result['pass'],'seconds',round(result['seconds'],1),flush=True)
    if not result['pass'] or not result['postrefinement_exclusion_pass']:return {'candidate_id':candidate,'status':'geometry_failed','audit':result}
    model=probe(dst,prefix,('A',),{i+1:r['human_pos'] for i,r in enumerate(target)})
    return {'candidate_id':candidate,'species':species,'status':'complete','minimum_pH_contrast':model['proton_model']['minimum_contrast_kcal']}

if __name__=='__main__':run(*sys.argv[1:])
