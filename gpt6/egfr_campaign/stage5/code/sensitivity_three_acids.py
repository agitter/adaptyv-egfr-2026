"""Expanded classical proton linkage: actual E106 contact plus D54/E103.
Corrects scope, not the earlier arithmetic. Three carboxylates plus two receptor
histidines give243 microstates. New outputs do not overwrite earlier81-state tests.
"""
from pathlib import Path
import sys,json,shutil
R=Path('/mnt/data/egfr_campaign');S=R/'stage5'
sys.path[:0]=[str(R/'stage3/code'),str(R/'stage2/code'),str(R/'code'),str(R/'stage2/runtime/python')]
import acid_histidine_ensemble as a
if __name__=='__main__':
 cid=sys.argv[1];root=S/'intermediate/three_acid_scope';(root/'intermediate/refined').mkdir(parents=True,exist_ok=True)
 source=S/'intermediate/refined'/f'{cid}_human6ARU_relaxed.pdb';dest=root/'intermediate/refined'/source.name;shutil.copy2(source,dest)
 d=json.loads((S/'intermediate/designs'/f'{cid}.json').read_text());pos=tuple(r['position'] for r in d['structure'] if r.get('parent_position') in [54,103,106]);assert len(pos)==3,pos
 a.S3=root;a.probe_acids(cid,pos)
