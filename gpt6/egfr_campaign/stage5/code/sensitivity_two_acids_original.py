"""Reuse classical carboxylate and terminal-cap sensitivity with correct indexing.
Deletion changes residue numbering: derive designed acid positions by parent map.
"""
from pathlib import Path
import sys,json
R=Path('/mnt/data/egfr_campaign');S=R/'stage5'
sys.path[:0]=[str(R/'stage3/code'),str(R/'stage2/code'),str(R/'code'),str(R/'stage2/runtime/python')]
def run(cid,mode):
 if mode=='acid':
  import acid_histidine_ensemble as a
  a.S3=S;d=json.loads((S/'intermediate/designs'/f'{cid}.json').read_text());pos=tuple(r['position'] for r in d['structure'] if r.get('parent_position') in [54,103]);assert len(pos)==2,pos
  a.probe_acids(cid,pos)
 elif mode=='cap':
  import cap_target as c
  c.S3=S;result=c.main(cid);print(cid,'capped minimum',result['proton_model']['minimum_contrast_kcal'],flush=True)
 else:raise ValueError(mode)
if __name__=='__main__':run(*sys.argv[1:])
