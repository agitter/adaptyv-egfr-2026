"""Use supplied antibody annotations as data, not a newly run annotator.
Keep original motif-only databases as auditable intermediate versions.
"""
import csv,gzip,json,re,hashlib,time
from classical_design import S,R
N=S/'intermediate/novelty';D=R/'inputs/transfer2/egfr_refinement_fetch'
csv.field_size_limit(32*1024*1024)
refs={r['sequence']:r['sources'] for r in map(json.loads,open(N/'cdrh3.jsonl'))};extra={};counts={'numbered_H_rows':0,'validated_annotated_cdr3':0,'anchor_or_length_rejected':0,'numbered_variable_domains':0}
for rownum,r in enumerate(csv.DictReader(gzip.open(D/'plabdab_unpaired.csv.gz','rt')),1):
 s=r.get('numbered','').upper();seq=s.replace('-','');src=f'plabdab_unpaired.csv.gz:row{rownum}:provided_numbered_variable_region'
 if len(seq)>=30 and set(seq)<=set('ACDEFGHIKLMNPQRSTVWYXBZ'):
  extra.setdefault(seq,src);counts['numbered_variable_domains']+=1
 if r['chain']!='H':continue
 counts['numbered_H_rows']+=1
 try:n=int(r['cdr_lengths'].split('_')[2])
 except Exception:continue
 if len(s)>=115 and s[103]=='C' and s[-11] in 'WF':
  h=s[104:-11].replace('-','')
  if len(h)==n and set(h)<=set('ACDEFGHIKLMNPQRSTVWYXBZ'):
   refs.setdefault(h,[])
   if len(refs[h])<3:refs[h].append(src+':validated_CDRH3')
   counts['validated_annotated_cdr3']+=1;continue
 counts['anchor_or_length_rejected']+=1
# Broaden conserved endpoint patterns; all windows are retained as conservative
# candidates rather than falsely claiming a definitive numbering assignment.
pat=re.compile(r'C([ACDEFGHIKLMNPQRSTVWYXBZ]{5,40}?)[WF][GA][ACDEFGHIKLMNPQRSTVWYXBZ][GA]')
for seq,src in extra.items():
 for i,c in enumerate(seq):
  if c=='C':
   m=pat.match(seq,i)
   if m:
    h=m.group(1);refs.setdefault(h,[])
    if len(refs[h])<3:refs[h].append(src+':broad_endpoint')
with open(N/'cdrh3_augmented.jsonl','w') as f:
 for seq,src in sorted(refs.items()):f.write(json.dumps({'sequence':seq,'sources':src})+'\n')
# Domain records can improve full-chain searches of antibodies with long constant
# regions. They are added as references, never as design seeds.
old={};key=None
for line in open(N/'antibodies.fasta'):
 if line.startswith('>'):key=line[1:].strip()
 elif key:old.setdefault(line.strip(),key)
with open(N/'antibodies_augmented.fasta','w') as f:
 for seq,key in old.items():f.write('>'+key+'\n'+seq+'\n')
 for seq,src in extra.items():
  if seq not in old:f.write('>numbered_'+hashlib.sha256(seq.encode()).hexdigest()[:20]+'\n'+seq+'\n')
counts.update({'cdr3_total':len(refs),'old_cdr3':sum(1 for _ in open(N/'cdrh3.jsonl')),'new_unique_domain_references':sum(seq not in old for seq in extra),'antibody_total':len(set(old)|set(extra))})
(S/'reference/novelty_augmentation.json').write_text(json.dumps(counts,indent=2));print(json.dumps(counts,indent=2))
