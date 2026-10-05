from pathlib import Path
import numpy as np,json,re,math
from Bio import Align,SeqIO
from Bio.PDB import MMCIFParser,PDBParser,ShrakeRupley
from Bio.SeqUtils import seq1
from scipy.spatial import cKDTree
R=Path('/mnt/data/egfr_campaign');O=R/'intermediate';D=R/'inputs/transfer1/legacy_fetch';P=R/'inputs/provided'
human=str(next(SeqIO.parse(O/'P00533.fasta','fasta')).seq); mouse=str(next(SeqIO.parse(O/'Q01279.fasta','fasta')).seq)
ali=Align.PairwiseAligner(mode='global',match_score=2,mismatch_score=-1,open_gap_score=-8,extend_gap_score=-0.5)
def mapping(ref,seq):
 a=ali.align(ref,seq)[0]; out={}
 for (a1,a2),(b1,b2) in zip(*a.aligned):
  for i,j in zip(range(a1,a2),range(b1,b2)):out[j]=i+1
 return out
hm=mapping(human,mouse); inv={v:k+1 for k,v in hm.items()}
np.save(O/'human_mouse_map.npy',np.array([[k,v] for k,v in inv.items()]))
print('ECR identity',sum(human[i-1]==mouse[inv[i]-1] for i in range(25,646) if i in inv)/621)
prepared={}
for name in ['human6ARU','mouseAF','1IVO','1NQL','1YY9']:
 a=json.loads((O/(name+'_chainA.json')).read_text()); seq=''.join(x['aa'] for x in a)
 mp=mapping(mouse if name=='mouseAF' else human,seq)
 for i,x in enumerate(a):
  x['uniprot']=mp[i]; x['human_pos']=(hm.get(mp[i]-1,hm.get(mp[i]-2,hm.get(mp[i]-3,0))+0.1)) if name=='mouseAF' else mp[i]
  x['conserved']=x['human_pos'] in inv and human[x['human_pos']-1]==mouse[inv[x['human_pos']]-1]
 a=[x for x in a if 25<=x['human_pos']<=645]
 prepared[name]=a
 print(name,'range',a[0]['uniprot'],a[-1]['uniprot'])
# Kabsch alignment of domain III core; coordinate arithmetic only, no prediction.
def fit(a,b):
 A=a-a.mean(0);B=b-b.mean(0);u,s,v=np.linalg.svd(A.T@B);m=np.eye(3);m[2,2]=np.linalg.det(u@v);rot=u@m@v;tr=b.mean(0)-a.mean(0)@rot
 return rot,tr,float(np.sqrt(np.mean(np.sum((a@rot+tr-b)**2,axis=1))))
base={x['human_pos']:x for x in prepared['human6ARU']};met={}
for name,a in prepared.items():
 ca=np.array([x['atoms']['CA'] for x in a if 334<=x['human_pos']<=504 and x['human_pos'] in base]);cb=np.array([base[x['human_pos']]['atoms']['CA'] for x in a if 334<=x['human_pos']<=504 and x['human_pos'] in base]); rot,tr,rms=fit(ca,cb)
 for x in a:x['atoms']={k:(np.array(v)@rot+tr).tolist() for k,v in x['atoms'].items()}
 (O/(name+'_aligned.json')).write_text(json.dumps(a))
 met[name]={'domainIII_alignment_rmsd_A':rms,'rotation':rot.tolist(),'translation':tr.tolist()};print(name,rms)
(R/'reference/target_alignment.json').write_text(json.dumps(met,indent=2))
# Coordinate-based solvent-accessible areas, 1973 rolling-probe sphere method.
parser=MMCIFParser(QUIET=True);st=parser.get_structure('human',str(P/'6ARU.cif'));sr=ShrakeRupley(n_points=200);sr.compute(st[0]['A'],level='R')
rr={r.id[1]:r for r in st[0]['A'] if r.id[0]==' '}
for a in prepared['human6ARU']:
 a['sasa']=float(rr[a['id']].sasa)
 a['mouse_aa']=mouse[inv[a['human_pos']]-1] if a['human_pos'] in inv else '-'
(O/'human6ARU_aligned.json').write_text(json.dumps(prepared['human6ARU']))
seqs={'human':human,'mouse':mouse};glyco={}
for sp,seq in seqs.items():
 glyco[sp]=[i+1 for i in range(24,643) if seq[i]=='N' and seq[i+1]!='P' and seq[i+2] in 'ST']
(R/'reference/glycosylation_sequons.json').write_text(json.dumps(glyco,indent=2))
print('Sequons',glyco)
print('Domain III exposed charged groups (UniProt positions):')
for a in prepared['human6ARU']:
 if 334<=a['human_pos']<=505 and a['aa'] in 'HDE' and a['sasa']>12:print(a['human_pos'],a['aa'],'mouse',a['mouse_aa'],'sasa',round(a['sasa'],1),'PDB',a['id'])
# Crystallographic glycan atoms, saved independently from protein.
gly=[]
for c in st[0]:
 for r in c:
  if r.resname in ['NAG','MAN','BMA','FUC','SIA','GAL','GLC']:
   for at in r:
    if at.element!='H':gly.append({'chain':c.id,'resid':r.id[1],'resname':r.resname,'atom':at.name,'xyz':at.coord.tolist()})
(O/'human6ARU_glycans.json').write_text(json.dumps(gly))
print('Glycan atoms',len(gly))
