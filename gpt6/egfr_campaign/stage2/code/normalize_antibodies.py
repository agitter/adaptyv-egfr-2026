"""Streaming sequence-only novelty reference. No model inputs are copied."""
from pathlib import Path
import csv,gzip,json,hashlib,re,ast,time
R=Path('/mnt/data/egfr_campaign');S=R/'stage2';D=R/'inputs/transfer2/egfr_refinement_fetch';out=S/'intermediate/novelty';out.mkdir(exist_ok=True)
csv.field_size_limit(32*1024*1024)
# Endpoint enumeration intentionally conservative; not an IMGT annotator.
pat=re.compile(r'C([ACDEFGHIKLMNPQRSTVWYXBZ]{5,40}?)[WF][GAS]QG')
seqs={};cdrs={};counts={}
def add(seq,src,annot=None):
 seq=re.sub(r'\s','',str(seq)).upper()
 if len(seq)<30 or not set(seq)<=set('ACDEFGHIKLMNPQRSTVWYXBZ'):return
 key=hashlib.sha256(seq.encode()).hexdigest()[:20]
 if key not in seqs:seqs[key]=[seq,[src]]
 elif len(seqs[key][1])<4:seqs[key][1].append(src)
 for i,c in enumerate(seq):
  if c!='C':continue
  m=pat.match(seq,i)
  if m:
   h=m.group(1);cdrs.setdefault(h,[])
   if len(cdrs[h])<3:cdrs[h].append(src)
 if annot:
  cdrs.setdefault(annot,[])
  if len(cdrs[annot])<3:cdrs[annot].append(src+':provided_CDRH3')
t0=time.time()
for fname in ['plabdab_paired.csv.gz','plabdab_unpaired.csv.gz','plabdab_nano.csv.gz','therasabdab.csv']:
 f=gzip.open(D/fname,'rt') if fname.endswith('.gz') else open(D/fname,encoding='utf-8-sig');n=0
 for row in csv.DictReader(f):
  n+=1
  if fname.startswith('plabdab_paired'):
   for k in ['heavy_sequence','light_sequence']:add(row.get(k,''),f'{fname}:row{n}:{k}')
  elif fname.startswith('plabdab_unpaired'):add(row.get('GBSeq_sequence',''),f'{fname}:row{n}')
  elif fname.startswith('plabdab_nano'):
   try:ann=ast.literal_eval(row.get('cdr_sequences','{}')).get('CDRH3')
   except Exception:ann=None
   add(row.get('sequence',''),f'{fname}:row{n}',ann)
  else:
   for k in ['HeavySequence','LightSequence','HeavySequence(ifbispec)','LightSequence(ifbispec)']:add(row.get(k,''),f'{fname}:row{n}:{k}')
  if n%100000==0:print(fname,n,'unique',len(seqs),'seconds',round(time.time()-t0,1),flush=True)
 f.close();counts[fname]=n;print('DONE',fname,n,len(seqs),len(cdrs),flush=True)
with open(out/'antibodies.fasta','w') as f,open(out/'antibody_sources.jsonl','w') as g:
 for key,(seq,src) in sorted(seqs.items()):f.write('>'+key+'\n'+seq+'\n');g.write(json.dumps({'id':key,'sources':src})+'\n')
# Include the earlier Swiss-Prot/PDB segments; they are database records only.
for line in open(R/'intermediate/known_cdrh3.jsonl'):
 r=json.loads(line);cdrs.setdefault(r['sequence'],r['sources'])
with open(out/'cdrh3.jsonl','w') as f:
 for seq,src in sorted(cdrs.items()):f.write(json.dumps({'sequence':seq,'sources':src})+'\n')
(S/'reference/novelty_reference_stats.json').write_text(json.dumps({'input_rows':counts,'unique_antibodies':len(seqs),'unique_cdr3':len(cdrs),'seconds':time.time()-t0,'annotation':'conservative endpoint enumeration plus provided nanobody CDR3 annotations; NOT ANARCI/IMGT'},indent=2))
