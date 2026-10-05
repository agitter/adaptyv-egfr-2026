"""Distinct-sequence catalog; rejected hypotheses are not counted as successes."""
from pathlib import Path
import json,hashlib,collections,re
R=Path('/mnt/data/egfr_campaign');S=R/'stage5'

def main():
 old={};counts={}
 for stage in ['stage2','stage3','stage4']:
  n=0
  for p in sorted((R/stage/'intermediate/designs').glob('*.json')):
   d=json.loads(p.read_text());seq=d.get('sequence')
   if seq:
    old.setdefault(seq,[]).append({'stage':stage,'candidate_id':d.get('candidate_id'),'path':str(p)});n+=1
  counts[stage]=n
 rows=[];seen={};duplicates=[];cd=[]
 for p in sorted((S/'intermediate/designs').glob('*.json')):
  d=json.loads(p.read_text())
  if 'sequence' not in d or d['candidate_id']=='B00000':continue
  seq=d['sequence'];rows.append({'candidate_id':d['candidate_id'],'sequence':seq,'length':len(seq),'molecule_class':d['molecule_class'],'branch':d['branch'],'cdr3':d['cdr_sequences'][2],'source':str(p),'source_sha256':hashlib.sha256(p.read_bytes()).hexdigest()});cd.append(d['cdr_sequences'][2])
  if seq in seen:duplicates.append([seen[seq],d['candidate_id']])
  seen[seq]=d['candidate_id']
  if seq in old:duplicates.append({'candidate_id':d['candidate_id'],'previous_records':old[seq]})
 ref=[]
 for p in sorted((S/'intermediate/backbones').glob('*.json')):
  d=json.loads(p.read_text())
  if 'stats' in d and 'attempts' in d['stats']:ref.append({'file':str(p),'attempts':d['stats']['attempts'],'accepted':len(d.get('accepted',[])),'statistics':d['stats']})
 result={'stage5_candidates':len(rows),'unique_stage5_sequences':len(seen),'unique_stage5_CDR3s':len(set(cd)),'length_range':[min(r['length'] for r in rows),max(r['length'] for r in rows)],'per_prefix':dict(collections.Counter(r['candidate_id'][0] for r in rows)),'historical_design_records':counts,'historical_distinct_sequences':len(old),'duplicates_within_or_before_stage5':duplicates,'backbone_attempts':sum(r['attempts'] for r in ref),'accepted_backbone_records':sum(r['accepted'] for r in ref),'backbone_searches':ref,'qualification':'Catalog of computational hypotheses including explicit failures; B00000 is an unchanged control, not a new candidate.'}
 (S/'reference/candidate_catalog.json').write_text(json.dumps(rows,indent=2));(S/'reference/catalog_summary.json').write_text(json.dumps(result,indent=2));assert not duplicates,duplicates;print({k:v for k,v in result.items() if k!='backbone_searches'},flush=True)
 # Every candidate is retained for scientific provenance, not implicitly recommended.
 out=S/'output';out.mkdir(exist_ok=True);(out/'stage5_276_NOT_FOR_SUBMISSION.fasta').write_text(''.join('>'+r['candidate_id']+' NOT_FOR_SUBMISSION_UNVALIDATED\n'+r['sequence']+'\n' for r in rows))
if __name__=='__main__':main()
