"""Unmodified historical finite cluster expansion, redirected to this round."""
from pathlib import Path
import sys
R=Path('/mnt/data/egfr_campaign');S=R/'stage7';sys.path[:0]=[str(R/x/'code') for x in ['stage6','stage5','stage3','stage2']]+[str(R/'code'),str(R/'stage2/runtime/python')]
import cluster_acids
from mm_refine import mm
mm.Platform.getPlatformByName('CPU').setPropertyDefaultValue('Threads','1')
cluster_acids.S=S
if __name__=='__main__':cluster_acids.main(sys.argv[1])
