"""Exact CDR3 edit screening against the already hashed, reconstructed references."""
from pathlib import Path
import json,hashlib,time
from rapidfuzz.distance import Levenshtein
from rapidfuzz import process
R=Path('/mnt/data/egfr_campaign');S=R/'stage5';O=S/'intermediate/novelty'

def main():
 start=time.time();p=O/'cdrh3_rebuilt.jsonl';sha=hashlib.sha256(p.read_bytes()).hexdigest();assert sha=='2268c5339550bb68acdcaed0ec7a72cc11f15bb0cf16ea6c68c5ea9d8417a436'
 refs={};by={}
 for line in p.open():
  r=json.loads(line);refs[r['sequence']]=r['sources'];by.setdefault(len(r['sequence']),[]).append(r['sequence'])
 assert len(refs)==203818
 cache={};rows=[]
 for p in sorted((S/'intermediate/designs').glob('H*.json')):
  d=json.loads(p.read_text());q=d['cdr_sequences'][2]
  if q not in cache:
   best=-1;winner=None;distance=None
   for length,seqs in by.items():
    if min(length,len(q))/max(length,len(q))<=best:continue
    hit=process.extractOne(q,seqs,scorer=Levenshtein.distance)
    if hit:
     ref,dist,index=hit;identity=1-dist/max(len(q),len(ref))
     if identity>best:best,winner,distance=identity,ref,dist
   cache[q]={'cdr3_query':q,'closest_sequence':winner,'closest_sources':refs[winner],'edit_distance':int(distance),'max_edit_identity':best,'below_70pct':best<.7}
  rows.append({'candidate_id':d['candidate_id'],**cache[q]})
 (O/'chemistry_cdr3_screen.json').write_text(json.dumps(rows,indent=2));summary={'candidate_count':len(rows),'unique_query_cdr3':len(cache),'reference_cdr3_count':len(refs),'reference_sha256':sha,'all_below_70pct':all(r['below_70pct'] for r in rows),'maximum_edit_identity':max(r['max_edit_identity'] for r in rows),'seconds':time.time()-start,'qualification':'Local historical edit search, not organizer antibody numbering or upload acceptance.'};(S/'reference/chemistry_novelty_summary.json').write_text(json.dumps(summary,indent=2));print(summary)
if __name__=='__main__':main()
