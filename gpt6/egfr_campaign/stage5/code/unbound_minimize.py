"""Same unbound1000-step-limit minimization for all E104 prototypes and control.
Local relaxation only; no claim of equilibrium folding or entropy."""
from pathlib import Path
import sys
R=Path('/mnt/data/egfr_campaign');S=R/'stage5'
sys.path[:0]=[str(R/'stage3/code'),str(R/'stage2/code'),str(R/'stage2/runtime/python'),str(R/'code')]
import free_binder
if __name__=='__main__':
 free_binder.S3=S;free_binder.main(sys.argv[1],0)
