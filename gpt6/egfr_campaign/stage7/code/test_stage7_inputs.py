"""Independent, deterministic input/lineage checks for the replacement experiment."""
from pathlib import Path
import json,hashlib,re,math,sys
import numpy as np
R=Path('/mnt/data/egfr_campaign');S=R/'stage7'
sys.path[:0]=[str(R/'stage2/code'),str(R/'code')]
from classical_design import LIB
from legacy_geometry import dihedral

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def main():
 checks=[]
 def check(name,condition,details=None):
  checks.append({'test':name,'pass':bool(condition),'details':details})
 parentp=R/'stage6/intermediate/designs/J00008.json';parent=json.loads(parentp.read_text())
 gen=json.loads((S/'reference/generation.json').read_text())['records'];ds=[json.loads((S/'intermediate/designs'/f"{r['candidate_id']}.json").read_text()) for r in gen if r['status']=='generated']
 check('23 records including one unchanged control',len(ds)==23)
 check('all records have distinct sequences',len({d['sequence'] for d in ds})==23)
 check('22 new sequences',sum(r.get('new_sequence',False) for r in gen)==22)
 check('M00000 is unchanged parent',ds[0]['candidate_id']=='M00000' and ds[0]['sequence']==parent['sequence'])
 for d in ds:
  cid=d['candidate_id'];seq=d['sequence'];rows=d['structure'];edit=d['replacement_experiment']['mutations'];ep={e['position'] for e in edit}
  check(cid+' standard residues and required length',len(seq)==125 and set(seq)<=set('ACDEFGHIKLMNPQRSTVWY'))
  check(cid+' sequence equals coordinate residue identities',seq==''.join(r['aa'] for r in rows))
  check(cid+' unchanged framework',all(a==b for a,b,row in zip(seq,parent['sequence'],rows) if not row['loop']))
  check(cid+' all edits are intended glycine replacements',all(e['from']=='G' and e['to'] in ['A','S'] and parent['sequence'][e['position']-1]=='G' and rows[e['position']-1]['loop'] for e in edit))
  check(cid+' no unrecorded edits',{i+1 for i,(a,b) in enumerate(zip(seq,parent['sequence'])) if a!=b}==ep)
  check(cid+' original framework CDR mask retained',[r['loop'] for r in rows]==[r['loop'] for r in parent['structure']])
  check(cid+' CDR intervals and strings consistent',d['cdr_intervals_zero_based']==parent['cdr_intervals_zero_based'] and d['cdr_sequences']==[seq[a:b] for a,b in d['cdr_intervals_zero_based']])
  check(cid+' source hashes match',sha(d['ancestry']['parent_file'])==d['ancestry']['parent_sha256'] and sha(d['ancestry']['backbone_file'])==d['ancestry']['backbone_sha256'])
  check(cid+' contact identities untouched',all(seq[i-1]==parent['sequence'][i-1] for i in [52,53,54,103,104,106]))
  check(cid+' no new cysteine or proline',all(a==b for a,b in zip(seq,parent['sequence']) if a in 'CP' or b in 'CP'))
  check(cid+' no glycosylation sequon',re.search(r'N[^P][ST]',seq) is None)
  check(cid+' inherited scores invalidated',d['metrics'].get('scores_valid') is False)
  check(cid+' all required heavy atoms present',all(set(['N','CA','C','O']+LIB[r['aa']]['names'])<=set(r['atoms']) for r in rows))
  for e in edit:
   i=e['position']-1;phi=math.degrees(dihedral(np.asarray(rows[i-1]['atoms']['C']),*[np.asarray(rows[i]['atoms'][n]) for n in ['N','CA','C']]))
   check(f'{cid} position {i+1} has negative phi',phi<0,phi)
 cdr=json.loads((S/'reference/novelty_summary.json').read_text())
 check('CDR reference SHA matches prior saved collection',cdr['reference_sha256']=='2268c5339550bb68acdcaed0ec7a72cc11f15bb0cf16ea6c68c5ea9d8417a436')
 check('all new CDR screens below local70percent',cdr['all_below_70pct'] and cdr['candidate_count']==23 and cdr['maximum_edit_identity']<=.5+1e-10)
 for db in ['swissprot','pdb','antibodies_augmented']:
  p=S/'reference'/f'panel_blast_{db}_summary.json';j=json.loads(p.read_text())
  print(db,'summary keys',list(j))
  check(db+' all candidate query IDs observed',j['candidate_queries']==23 and j['observed_queries_including_control']==24 and j['all_query_ids_verified'])
  check(db+' exact reference positive control',j['positive_control_exact_match'])
  cmd=j['command'];q=Path(cmd[cmd.index('-query')+1]);o=Path(cmd[cmd.index('-out')+1])
  check(db+' query SHA matches observed record',sha(q)==j['query_sha256'])
  check(db+' result SHA matches observed record',sha(o)==j['result_sha256'])
 out={'checks':checks,'passed':sum(c['pass'] for c in checks),'count':len(checks),'all_pass':all(c['pass'] for c in checks),'qualification':'Input and implementation checks, not biological validation.'}
 (S/'reference/input_validation.json').write_text(json.dumps(out,indent=2))
 print('INPUT_CHECKS',out['passed'],'/',out['count'],flush=True)
 if not out['all_pass']:
  print(json.dumps([c for c in checks if not c['pass']],indent=2));raise AssertionError('Input audit failure')
if __name__=='__main__':main()
