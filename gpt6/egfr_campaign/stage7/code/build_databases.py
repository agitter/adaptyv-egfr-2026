from pathlib import Path
import sys,importlib.util
R=Path('/mnt/data/egfr_campaign');S=R/'stage7'
spec=importlib.util.spec_from_file_location('build_old',R/'stage5/code/build_databases.py');b=importlib.util.module_from_spec(spec);spec.loader.exec_module(b)
b.S=S;b.OUT=S/'intermediate/blast';b.BIN=S/'runtime/bin';b.OUT.mkdir(parents=True,exist_ok=True);b.BIN.mkdir(parents=True,exist_ok=True);b.main()
