from pathlib import Path
import sys
R=Path('/mnt/data/egfr_campaign');S=R/'stage6';sys.path.insert(0,str(R/'stage5/code'))
import build_databases as b
b.S=S;b.OUT=S/'intermediate/blast';b.BIN=S/'runtime/bin'
b.OUT.mkdir(parents=True,exist_ok=True);b.BIN.mkdir(parents=True,exist_ok=True)
b.main()
