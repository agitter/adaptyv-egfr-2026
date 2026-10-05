from pathlib import Path
import json,numpy as np,re,hashlib,time
from Bio import SeqIO
from legacy_geometry import atom_frame,dihedral
R=Path('/mnt/data/egfr_campaign');O=R/'intermediate'
expected={'A':'CB','C':'CB SG','D':'CB CG OD1 OD2','E':'CB CG CD OE1 OE2','F':'CB CG CD1 CD2 CE1 CE2 CZ','G':'','H':'CB CG ND1 CD2 CE1 NE2','I':'CB CG1 CG2 CD1','K':'CB CG CD CE NZ','L':'CB CG CD1 CD2','M':'CB CG SD CE','N':'CB CG OD1 ND2','P':'CB CG CD','Q':'CB CG CD OE1 NE2','R':'CB CG CD NE CZ NH1 NH2','S':'CB OG','T':'CB OG1 CG2','V':'CB CG1 CG2','W':'CB CG CD1 CD2 NE1 CE2 CE3 CZ2 CZ3 CH2','Y':'CB CG CD1 CD2 CE1 CE2 CZ OH'}
pool={aa:[] for aa in expected}
for name in ['human6ARU','1IVO','1NQL','1YY9','2A3D','1FSD','1L2Y']:
 rows=json.loads((O/(name+'_chainA.json')).read_text())
 for r in rows:
  aa=r['aa'];a=r['atoms'];names=expected[aa].split()
  if aa=='G' or not all(k in a for k in ['N','CA','C']+names):continue
  frame=atom_frame(*[np.array(a[k],float) for k in ['N','CA','C']]);xyz=(np.array([a[k] for k in names])-a['CA'])@frame.T
  pool[aa].append({'xyz':xyz,'source':name+':A:'+str(r['id'])})
lib={}
for aa,names in expected.items():
 if aa=='G':lib[aa]={'names':[],'xyz':[[]],'sources':['empty glycine side chain']};continue
 choices=[]
 for item in pool[aa]:
  if not choices or min(np.sqrt(np.mean((item['xyz']-q['xyz'])**2)) for q in choices)>.6:choices.append(item)
  if len(choices)>=18:break
 if not choices:raise ValueError('Missing amino acid '+aa)
 lib[aa]={'names':names.split(),'xyz':[q['xyz'].tolist() for q in choices],'sources':[q['source'] for q in choices]}
 print(aa,len(pool[aa]),len(choices))
(O/'empirical_rotamers.json').write_text(json.dumps(lib))
# Extract conservative heavy-chain CDR3 candidates by conserved Cys/FR4 anchors.
# This is a motif-based coverage screen, not an IMGT annotation claim.
pat=re.compile(r'C([ACDEFGHIKLMNPQRSTVWY]{5,35})W[GAS]QG')
unique={};counts={}
for name in ['swissprot','pdb']:
 n=0;matches=0
 for rec in SeqIO.parse(O/(name+'.fasta'),'fasta'):
  n+=1;s=str(rec.seq)
  for m in pat.finditer(s):
   c=m.group(1)
   if c not in unique:unique[c]=[name+':'+rec.id]
   elif len(unique[c])<3:unique[c].append(name+':'+rec.id)
   matches+=1
 counts[name]={'sequences':n,'motif_matches':matches};print(name,counts[name],flush=True)
with open(O/'known_cdrh3.jsonl','w') as f:
 for seq,src in sorted(unique.items()):f.write(json.dumps({'sequence':seq,'sources':src})+'\n')
counts['unique_cdrh3']=len(unique)
(R/'reference/cdr_database_stats.json').write_text(json.dumps(counts,indent=2))
print('Unique CDR3 candidates',len(unique))
