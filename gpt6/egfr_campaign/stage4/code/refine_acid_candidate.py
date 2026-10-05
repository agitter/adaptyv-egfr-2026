"""Second, explicitly protonated structural hypothesis using historical MM."""
from pathlib import Path
import sys
R=Path('/mnt/data/egfr_campaign');S=R/'stage4'
sys.path[:0]=[str(R/'stage3/code'),str(R/'stage2/code'),str(R/'code'),str(R/'stage2/runtime/python')]
import refine_state
if __name__=='__main__':
 refine_state.S3=S
 print(refine_state.run(sys.argv[1],sys.argv[2] if len(sys.argv)>2 else 'human6ARU',600),flush=True)
