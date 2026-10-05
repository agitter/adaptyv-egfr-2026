"""Launch a bounded synchronous-campaign task with preserved configuration."""
from pathlib import Path
import os,subprocess,json,sys,datetime
S=Path('/mnt/data/egfr_campaign/stage4');R=S.parent
cfg={'candidate':'S00022','seed':198003,'production_ps':20.,'warmup_ps':2.,'threads':2,'rationale':'Disulfide did not outperform equal-refinement parent in extended trajectory. S00022 has no extra cysteines and passed both mouse minimized geometry and short target-overlay audit. Apply identical extended protocol, without choosing a seed based on results.','created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
(S/'reference/final_proline_predeclaration.json').write_text(json.dumps(cfg,indent=2))
env=os.environ.copy();env.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',OPENMM_CPU_THREADS='1');env['PYTHONPATH']=':'.join(str(p) for p in [S/'code',R/'stage3/code',R/'stage2/code',R/'code',R/'stage2/runtime/python'])
cmd=[sys.executable,str(S/'code/final_proline_comparison.py')]
with open(S/'logs/final_proline_comparison.log','w') as f:p=subprocess.Popen(cmd,env=env,stdout=f,stderr=subprocess.STDOUT,start_new_session=True)
(S/'reference/final_proline_pid.json').write_text(json.dumps({'pid':p.pid,'command':cmd,'must_finish_before_checkpoint':True},indent=2));print(p.pid)
