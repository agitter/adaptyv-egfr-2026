"""Matched classical interface accessibility (Shrake-Rupley1973), not affinity."""
from pathlib import Path
import sys
R=Path('/mnt/data/egfr_campaign');S=R/'stage4'
sys.path[:0]=[str(R/'stage3/code'),str(R/'stage2/code'),str(R/'stage2/runtime/python')]
import interface_surface
if __name__=='__main__':
 for root,cid in [(S,'S00022'),(S,'S00028'),(S/'intermediate/matched_controls','C00003')]:
  interface_surface.S=root;interface_surface.main(cid,192)
