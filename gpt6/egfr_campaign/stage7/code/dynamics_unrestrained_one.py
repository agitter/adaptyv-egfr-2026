"""Classical matched thermal test with NO added chirality-restraint potential.
Standard Amber bonded/improper terms and HBond constraints remain. The inherited
warmup positional restraint ramps to zero before production. Job-isolated input
copies prevent concurrent native-control preparation from sharing output files.
"""
from pathlib import Path
import sys,json,hashlib,shutil,datetime
R=Path('/mnt/data/egfr_campaign');S=R/'stage7'
sys.path[:0]=[str(R/'stage4/code'),str(R/'stage3/code'),str(R/'stage2/code'),str(R/'code'),str(R/'stage2/runtime/python')]
import dynamics
from mm_refine import mm
mm.Platform.getPlatformByName('CPU').setPropertyDefaultValue('Threads','1')
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
cid=sys.argv[1];seed=int(sys.argv[2]);job=S/'intermediate/thermal_unrestrained'/f'{cid}_s{seed}'
for sub in ['intermediate/designs','intermediate/refined','intermediate/dynamics','reference']:(job/sub).mkdir(parents=True,exist_ok=True)
manifest=[]
if cid!='NATIVE_3EAK':
 for sub,name in [('designs',cid+'.json'),('refined',cid+'_human6ARU_relaxed.pdb')]:
  src=S/'intermediate'/sub/name;dst=job/'intermediate'/sub/name
  if dst.exists() and sha(dst)!=sha(src):raise ValueError('Different input already exists')
  shutil.copyfile(src,dst);manifest.append({'source':str(src),'copy':str(dst),'sha256':sha(src)})
else:
 src=R/'inputs/transfer1/legacy_fetch/3EAK.pdb';manifest.append({'source':str(src),'sha256':sha(src),'role':'native fold control only'})
(job/'reference/input_copies.json').write_text(json.dumps(manifest,indent=2))
(job/'reference/protocol.json').write_text(json.dumps({'candidate_id':cid,'seed':seed,'started_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'engine_file':str(R/'stage4/code/dynamics.py'),'engine_sha256':sha(R/'stage4/code/dynamics.py'),'extra_chirality_restraints':False,'standard_forcefield_impropers_retained':True,'production_ps':20.0,'warmup_ps':2.0,'threads':1,'interpretation':'Local thermal stress only; not an equilibrium stability or folding calculation.'},indent=2))
dynamics.S=job
dynamics.run(cid,seed,20.,2.,threads=1,chirality=False)
result=job/'intermediate/dynamics'/f'{cid}_v2_s{seed}.json';d=json.loads(result.read_text())
assert d['chirality_restraints'] is False and d['position_restraints_production'] is False
print('NO_ADDED_CHIRALITY_RESTRAINTS_CONFIRMED',cid,seed,flush=True)
