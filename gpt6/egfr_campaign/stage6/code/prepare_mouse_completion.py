"""Prepare mouse follow-ups of six already generated, human-geometry-passing designs.
Use each saved relaxed human binder as the starting geometry. Sequences are not
changed and are not counted as new designs.
"""
from pathlib import Path
import json,sys,hashlib
import numpy as np
R=Path('/mnt/data/egfr_campaign');S=R/'stage6';sys.path[:0]=[str(R/'stage2/code'),str(R/'stage2/runtime/python'),str(R/'code')]
from mm_refine import app,u
from Bio.SeqUtils import seq1
keys=['stage4:S00048','stage5:H00018','stage5:T3080600','stage5:T3070900','stage5:H00007','stage5:H00016'];models=json.loads((S/'reference/uniform_model_summary.json').read_text());selected=[];jobs=[]
for i,key in enumerate(keys,1):
 m=next(r for r in models if r['key']==key and r['species']=='human' and not r['acid_refined'] and r['combined_geometry_pass']);dp=Path(m['design_source']);d=json.loads(dp.read_text());p=app.PDBFile(m['coordinate_source']);xyz=np.asarray(p.positions.value_in_unit(u.angstrom));res=[r for r in p.topology.residues() if r.chain.id=='B'];sequence=''.join(seq1(r.name,custom_map={'HID':'H','HIE':'H','HIP':'H','CYX':'C','ASH':'D','GLH':'E'}) for r in res);assert sequence==d['sequence']
 for row,r in zip(d['structure'],res):row['atoms']={a.name:xyz[a.index].tolist() for a in r.atoms() if a.element!=app.element.hydrogen and a.name!='OXT'}
 cid=f'L{i:05d}';d.update(candidate_id=cid,campaign_stage=6,reassessment_source={'candidate':key,'file':str(dp),'sha256':hashlib.sha256(dp.read_bytes()).hexdigest(),'sequence_changed':False,'relaxed_human_starting_coordinates':m['coordinate_source'],'coordinate_sha256':hashlib.sha256(Path(m['coordinate_source']).read_bytes()).hexdigest()},metrics={'scores_valid':False,'note':'Previously evaluated sequence, fresh mouse follow-up; no inherited score substitution.'});(S/'intermediate/designs'/f'{cid}.json').write_text(json.dumps(d));selected.append({'candidate_id':cid,'source_key':key,'sequence':sequence});jobs.append({'id':cid+'_mouseAF','command':[sys.executable,str(S/'code/refine_one.py'),cid,'mouseAF'],'timeout':600})
plan={'name':'mouse_completion','concurrency':2,'selection':'Six existing candidates with passing detailed human geometry and positive human pH contrast; before these mouse results. No sequence generation.','candidates':selected,'jobs':jobs};(S/'reference/mouse_completion_plan.json').write_text(json.dumps(plan,indent=2));print('Prepared',len(jobs),'mouse follow-ups; not launched')
