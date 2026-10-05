"""Declared extension after all shorter contact-preserving returns failed.
Retains failed length4/5/6 attempts; explores length7 and unchanged-length8.
"""
from pathlib import Path
import sys,json,time
S=Path('/mnt/data/egfr_campaign/stage5');sys.path.insert(0,str(S/'code'))
import build_preserved_tail as b
import design_preserved_tail as d
if __name__=='__main__':
 for n in [7,8]:
  result=b.run('T3',n,1800,16)
  if result['stats']['diverse']:d.run(f'T3_n{n:02d}',3)
