"""Preserved velocity-Verlet/RATTLE and MT19937 Andersen-style protocol.
20 ps production, 2 ps ramp. This is a local stress test, not folding validation.
"""
from pathlib import Path
import sys,json,time
R=Path('/mnt/data/egfr_campaign');S=R/'stage7'
sys.path[:0]=[str(R/x/'code') for x in ['stage4','stage3','stage2']]+[str(R/'code'),str(R/'stage2/runtime/python')]
import dynamics
from mm_refine import mm
mm.Platform.getPlatformByName('CPU').setPropertyDefaultValue('Threads','1')
dynamics.S=S
if __name__=='__main__':dynamics.run(sys.argv[1],int(sys.argv[2]),20.,2.,threads=1,chirality=True)
