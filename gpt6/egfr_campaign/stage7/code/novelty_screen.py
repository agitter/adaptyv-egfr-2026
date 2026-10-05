"""Streaming reconstruction of supplied CDR3 references, then exact edit search.
No learned numbering or sequence model. Supplied modern annotations are data.
"""
from pathlib import Path
import csv,gzip,json,re,ast,time,sys,hashlib
from rapidfuzz.distance import Levenshtein
from rapidfuzz import process
R=Path('/mnt/data/egfr_campaign');S=R/'stage7';D=R/'inputs/transfer2/egfr_refinement_fetch';OUT=S/'intermediate/novelty';OUT.mkdir(parents=True,exist_ok=True)
csv.field_size_limit(32*1024*1024);AA=set('ACDEFGHIKLMNPQRSTVWYXBZ');strict=re.compile(r'C([ACDEFGHIKLMNPQRSTVWYXBZ]{5,40}?)[WF][GAS]QG');broad=re.compile(r'C([ACDEFGHIKLMNPQRSTVWYXBZ]{5,40}?)[WF][GA][ACDEFGHIKLMNPQRSTVWYXBZ][GA]')

def main():
 start=time.time();refs={};counts={}
 def add_cdr(seq,source):
  refs.setdefault(seq,[])
  if len(refs[seq])<3 and source not in refs[seq]:refs[seq].append(source)
 def scan(seq,source,pattern=strict):
  seq=re.sub(r'\s','',str(seq)).upper()
  if len(seq)<30 or not set(seq)<=AA:return
  for i,c in enumerate(seq):
   if c=='C':
    m=pattern.match(seq,i)
    if m:add_cdr(m.group(1),source)
 for fname in ['plabdab_paired.csv.gz','plabdab_unpaired.csv.gz','plabdab_nano.csv.gz','therasabdab.csv']:
  stream=gzip.open(D/fname,'rt') if fname.endswith('.gz') else open(D/fname,encoding='utf-8-sig');n=0
  for n,row in enumerate(csv.DictReader(stream),1):
   src=f'{fname}:row{n}'
   if fname=='plabdab_paired.csv.gz':
    for k in ['heavy_sequence','light_sequence']:scan(row.get(k,''),src+':'+k)
   elif fname=='plabdab_unpaired.csv.gz':
    scan(row.get('GBSeq_sequence',''),src)
    numbered=row.get('numbered','').upper();seq=numbered.replace('-','');scan(seq,src+':numbered_broad_endpoint',broad)
    if row.get('chain')=='H':
     try:length=int(row['cdr_lengths'].split('_')[2])
     except Exception:length=-1
     if len(numbered)>=115 and numbered[103]=='C' and numbered[-11] in 'WF':
      h=numbered[104:-11].replace('-','')
      if len(h)==length and set(h)<=AA:add_cdr(h,src+':provided_numbered_validated_CDRH3')
   elif fname=='plabdab_nano.csv.gz':
    sequence=re.sub(r'\s','',row.get('sequence','')).upper();scan(sequence,src)
    if len(sequence)>=30 and set(sequence)<=AA:
     try:ann=ast.literal_eval(row.get('cdr_sequences','{}')).get('CDRH3')
     except Exception:ann=None
     if ann:add_cdr(ann,src+':provided_CDRH3')
   else:
    for k in ['HeavySequence','LightSequence','HeavySequence(ifbispec)','LightSequence(ifbispec)']:scan(row.get(k,''),src+':'+k)
  stream.close();counts[fname]=n;print('SCANNED',fname,n,'CDR3',len(refs),flush=True)
 for line in open(R/'intermediate/known_cdrh3.jsonl'):
  r=json.loads(line)
  for src in r['sources'][:3]:add_cdr(r['sequence'],src)
 cache_path=OUT/'cdrh3_rebuilt.jsonl'
 with cache_path.open('w') as f:
  for seq,src in sorted(refs.items()):f.write(json.dumps({'sequence':seq,'sources':src})+'\n')
 by={}
 for seq in refs:by.setdefault(len(seq),[]).append(seq)
 cache={};records=[]
 for p in sorted((S/'intermediate/designs').glob('M*.json')):
  d=json.loads(p.read_text());q=d['cdr_sequences'][2]
  if q not in cache:
   best=-1;bestseq=None;distance=None
   for length,seqs in by.items():
    if min(length,len(q))/max(length,len(q))<=best:continue
    hit=process.extractOne(q,seqs,scorer=Levenshtein.distance)
    if hit:
     ref,dist,idx=hit;identity=1-dist/max(len(q),len(ref))
     if identity>best:best,bestseq,distance=identity,ref,dist
   cache[q]={'cdr3_query':q,'closest_sequence':bestseq,'closest_sources':refs[bestseq],'edit_distance':int(distance),'max_edit_identity':best,'below_70pct':best<.7}
  records.append({'candidate_id':d['candidate_id'],**cache[q]})
 (OUT/'stage7_cdr3_screen.json').write_text(json.dumps(records,indent=2))
 result={'input_rows':counts,'reference_cdr3_count':len(refs),'prior_expected_count':203818,'count_matches_prior':len(refs)==203818,'reference_sha256':hashlib.sha256(cache_path.read_bytes()).hexdigest(),'candidate_count':len(records),'unique_query_cdr3':len(cache),'all_below_70pct':all(r['below_70pct'] for r in records),'maximum_edit_identity':max(r['max_edit_identity'] for r in records),'seconds':time.time()-start,'scope':'CDR3 screening only; separately recorded BLAST searches provide whole-chain coverage when complete. Not official organizer annotation.'};(S/'reference/novelty_summary.json').write_text(json.dumps(result,indent=2));print(result,flush=True)
if __name__=='__main__':main()
