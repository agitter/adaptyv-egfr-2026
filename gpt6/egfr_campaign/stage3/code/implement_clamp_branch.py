from pathlib import Path
import hashlib,json
s=Path('/mnt/data/egfr_campaign/stage3');p=s/'code/gated_design.py';old=p.read_text()
(s/'reference/gated_design_preclamp.py.txt').write_text(old)
a=old.replace('def candidate_table(rec,allow_histidine=False,maxoptions=36):','def candidate_table(rec,allow_histidine=False,maxoptions=36,clamp_branch=False):')
a=a.replace("cache.append({'aa':aa", "if clamp_branch:\n                    from clamp_geometry import option_features\n                    pre-=5.0*float(option_features(aa,atoms).sum())\n                cache.append({'aa':aa")
a=a.replace("if counts.get(c['aa'],0)>=(4 if c['aa'] in 'DEYH' else 2):continue", "if counts.get(c['aa'],0)>=((9 if c['aa'] in 'DE' else 4) if clamp_branch else (4 if c['aa'] in 'DEYH' else 2)):continue")
assert a!=old;p.write_text(a)
(s/'reference/clamp_patch.json').write_text(json.dumps({'old_sha256':hashlib.sha256(old.encode()).hexdigest(),'new_sha256':hashlib.sha256(a.encode()).hexdigest(),'purpose':'New optional two-nitrogen branch; original default calculations unchanged.'},indent=2))
