"""Joint target/binder refinement, preserving the receptor during interface repair."""
from mm_binding_probe import *

def run(candidate,species='human6ARU',iterations=500):
    out=S/'intermediate/mm_binding';out.mkdir(exist_ok=True)
    d=json.loads((S/'intermediate/designs'/(candidate+'.json')).read_text());rows=d['structure'];target=json.loads((O/(species+'_aligned.json')).read_text());target=complete_target([r for r in target if 334<=r['human_pos']<=505]);assert len(target)==172
    tp=out/(candidate+'_'+species+'_target.pdb');bp=out/(candidate+'_binder_initial.pdb');write_pdb(target,tp,'A');write_pdb(rows,bp,'B')
    a=app.PDBFile(str(tp));b=app.PDBFile(str(bp));mod=app.Modeller(a.topology,a.positions);mod.add(b.topology,b.positions)
    src=out/(candidate+'_'+species+'_complex_initial.pdb');dst=out/(candidate+'_'+species+'_complex_relaxed.pdb')
    with open(src,'w') as f:app.PDBFile.writeFile(mod.topology,mod.positions,f)
    gly=json.loads((S/'intermediate/target_glycans.json').read_text()) if species=='human6ARU' else []
    mask=[False]*len(target)+[r['loop'] for r in rows]
    result=minimize(src,dst,mask,iterations,[g['xyz'] for g in gly]);result.update({'candidate_id':candidate,'species':species,'target_crop':[334,505],'joint_refinement':True,'glycan_obstacle_count':len(gly)})
    (out/(candidate+'_'+species+'_complex_audit.json')).write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
if __name__=='__main__':run(*sys.argv[1:])
