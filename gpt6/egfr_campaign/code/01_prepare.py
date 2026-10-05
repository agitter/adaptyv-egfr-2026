from pathlib import Path
import json,gzip,sys,subprocess,hashlib,platform,os
import numpy as np
from Bio.PDB import MMCIFParser,PDBParser
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.SeqUtils import seq1
from Bio import Align,SeqIO
R=Path('/mnt/data/egfr_campaign'); P=R/'inputs/provided'; D=R/'inputs/transfer1/legacy_fetch'; O=R/'intermediate'
mp=MMCIFParser(QUIET=True); pp=PDBParser(QUIET=True)
records={}
for name,fn,parser in [('human6ARU',P/'6ARU.cif',mp),('mouseAF',P/'AF-Q01279-F1-model_v6.cif',mp)]+[(s,D/(s+'.pdb'),pp) for s in ['1IVO','1NQL','1YY9','2A3D','1FSD','1L2Y','3DWT','3EAK','3EBA']]:
 st=parser.get_structure(name,str(fn)); chains=[]
 for c in st[0]:
  rr=[r for r in c if r.id[0]==' ' and 'CA' in r and seq1(r.resname)!='X']
  if not rr:continue
  seq=''.join(seq1(r.resname) for r in rr)
  chains.append({'chain':c.id,'nres':len(rr),'seq':seq,'first':str(rr[0].id),'last':str(rr[-1].id)})
  if c.id=='A':
   obj=[]
   for r in rr:
    obj.append({'id':r.id[1],'icode':r.id[2].strip(),'aa':seq1(r.resname),'atoms':{a.name:[float(x) for x in a.coord] for a in r if a.element!='H'},'bfactor':float(r['CA'].bfactor)})
   (O/(name+'_chainA.json')).write_text(json.dumps(obj))
 records[name]=chains
 print(name,[(x['chain'],x['nres'],x['first'],x['last']) for x in chains])
 if name in ['3EAK','3DWT','3EBA']:
  print(chains[0]['seq'])
  print(' '.join(str(r.id[1])+r.id[2].strip()+':'+seq1(r.resname) for r in list(st[0]['A']) if r.id[0]==' '))
(R/'reference/structure_inventory.json').write_text(json.dumps(records,indent=2))
with gzip.open(D/'uniprot_sprot.fasta.gz','rt') as f:
 for rec in SeqIO.parse(f,'fasta'):
  if rec.id.split('|')[1] in ['P00533','Q01279']:
   (O/(rec.id.split('|')[1]+'.fasta')).write_text('>'+rec.id+'\n'+str(rec.seq)+'\n')
   print(rec.id,len(rec.seq))
# Decompress with streaming for BLAST. No sequence contents held collectively in RAM.
import shutil
for src,dst in [('uniprot_sprot.fasta.gz','swissprot.fasta'),('pdb_seqres.txt.gz','pdb.fasta')]:
 with gzip.open(D/src,'rb') as a, open(O/dst,'wb') as b:shutil.copyfileobj(a,b)
 print(dst,(O/dst).stat().st_size)
