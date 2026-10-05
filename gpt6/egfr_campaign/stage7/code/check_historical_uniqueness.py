"""Reconcile new sequences against every retained prior candidate record.
This is within-campaign exact identity, not the organizer novelty calculation.
"""
from pathlib import Path
import json,hashlib
R=Path('/mnt/data/egfr_campaign');S=R/'stage7'
def main():
 history={};sources=[];skipped=[]
 for stage in ['stage2','stage3','stage4','stage5','stage6']:
  for p in sorted((R/stage/'intermediate/designs').glob('*.json')):
   raw=p.read_bytes();d=json.loads(raw)
   if not isinstance(d,dict) or 'sequence' not in d:skipped.append(str(p));continue
   seq=d['sequence'];history.setdefault(seq,[]).append(str(p));sources.append({'path':str(p),'sha256':hashlib.sha256(raw).hexdigest(),'sequence_sha256':hashlib.sha256(seq.encode()).hexdigest()})
 rows=[]
 for p in sorted((S/'intermediate/designs').glob('M*.json')):
  d=json.loads(p.read_text());hits=history.get(d['sequence'],[]);rows.append({'candidate_id':d['candidate_id'],'same_sequence_prior_records':hits,'is_new_to_retained_history':not hits})
 assert len(rows)==23 and sum(r['is_new_to_retained_history'] for r in rows)==22
 assert next(r for r in rows if r['candidate_id']=='M00000')['same_sequence_prior_records']
 out={'historical_candidate_records':len(sources),'historical_unique_sequences':len(history),'new_unique_sequences':22,'records':rows,'source_manifest':sources,'non_candidate_files_skipped':skipped,'qualification':'Exact identity against retained campaign designs only. Not a global novelty or official acceptance result.'}
 (S/'reference/historical_sequence_uniqueness.json').write_text(json.dumps(out,indent=2));print('Historical records',len(sources),'unique',len(history),'new sequences',22)
if __name__=='__main__':main()
