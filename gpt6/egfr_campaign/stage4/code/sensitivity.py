"""Classical charge-state and cropped-terminus sensitivity checks."""
from pathlib import Path
import sys,time
R=Path('/mnt/data/egfr_campaign');S=R/'stage4'
sys.path[:0]=[str(R/'stage3/code'),str(R/'stage2/code'),str(R/'code'),str(R/'stage2/runtime/python')]
def main(cid,mode):
 root=S/'intermediate/matched_controls' if cid.startswith('MATCHED_') else S;cid=cid.removeprefix('MATCHED_');p=root/'intermediate/refined'/(cid+'_human6ARU_proton.json');start=time.time()
 while not p.exists():
  if time.time()-start>600:raise TimeoutError('Source refinement unavailable')
  time.sleep(2)
 if mode=='acid':
  import acid_histidine_ensemble as a
  a.S3=root;a.probe_acids(cid,(54,103))
 elif mode=='cap':
  import cap_target as a
  a.S3=root;r=a.main(cid);print(cid,r['proton_model']['minimum_contrast_kcal'],flush=True)
 else:raise ValueError(mode)
if __name__=='__main__':main(*sys.argv[1:])
