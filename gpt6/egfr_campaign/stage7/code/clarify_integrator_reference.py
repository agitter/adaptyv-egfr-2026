"""Correct a bibliographic date in the report template, not any computation."""
from pathlib import Path
import json,hashlib,shutil
S=Path('/mnt/data/egfr_campaign/stage7')
p=S/'code/write_revision_report.py'
old=p.read_text()
needle='velocity-Verlet/RATTLE(1967/1983)'
assert old.count(needle)==1
backup=S/'reference/write_revision_report_before_integrator_date_clarification.py'
assert not backup.exists()
shutil.copy2(p,backup)
new=old.replace(needle,'velocity-Verlet/RATTLE(1982/1983)')
p.write_text(new)
r={'scope':'Report bibliographic date only. Verlet 1967 concerns the original position formulation; the velocity-Verlet reference is Swope et al.1982. No scientific code, model, trajectory or selection threshold was changed.','old_report_template_sha256':hashlib.sha256(old.encode()).hexdigest(),'new_report_template_sha256':hashlib.sha256(new.encode()).hexdigest(),'velocity_Verlet_reference':{'authors':'Swope WC, Andersen HC, Berens PH, Wilson KR','year':1982,'journal':'Journal of Chemical Physics76:637-649','doi':'10.1063/1.442716','title':'A computer simulation method for the calculation of equilibrium constants for the formation of physical clusters of molecules: Application to small water clusters','publisher_access':'DOI resolved to publisher abstract URL, but web retrieval failed; no new PDF or full-text analysis was performed.'},'RATTLE_reference':{'author':'Andersen HC','year':1983,'journal':'Journal of Computational Physics52:24-34','doi':'10.1016/0021-9991(83)90014-1','title':'Rattle: A velocity version of the shake algorithm for molecular dynamics calculations','publisher_abstract_index':'https://www.sciencedirect.com/science/article/pii/0021999183900141'}}
(S/'reference/integrator_bibliographic_clarification.json').write_text(json.dumps(r,indent=2))
print('Integrator bibliographic date clarified; no numerical changes.')
