"""Run the audited pre-2011 thermal protocol on stage5 inputs without changing it.
Two seeds and native/unchanged controls are predeclared. Dynamics are short
stress tests, not a prediction of folding free energy, kinetics or affinity.
"""
from pathlib import Path
import sys,json
R=Path('/mnt/data/egfr_campaign');S=R/'stage5'
sys.path[:0]=[str(R/'stage4/code'),str(R/'stage3/code'),str(R/'stage2/code'),str(R/'code'),str(R/'stage2/runtime/python')]
import dynamics
if __name__=='__main__':
 dynamics.S=S
 dynamics.run(sys.argv[1],int(sys.argv[2]),20.,2.,threads=2,chirality=True)
