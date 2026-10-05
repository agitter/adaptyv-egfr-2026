"""Build antigen-binding loops from random non-binder torsion fragments.
No original antibody CDR sequences or CDR coordinates initialize the loops.
Only generic framework atoms and loop attachment anchors are retained.
"""
from pathlib import Path
import numpy as np,json,sys,time
from legacy_geometry import *
R=Path('/mnt/data/egfr_campaign');O=R/'intermediate'; out=R/'stage2/intermediate/loops';out.mkdir(exist_ok=True)
source=json.loads((O/'3EAK_chainA.json').read_text());byid={x['id']:x for x in source}
regions={'H1':(26,38),'H2':(53,67),'H3':(100,117)}
mask=set(i for lo,hi in regions.values() for i in range(lo,hi+1))
framework=[r for r in source if r['id'] not in mask]
assert json.loads((O/'framework_only.json').read_text()) == framework
(R/'stage2/reference/framework_definition.json').write_text(json.dumps({'source':'3EAK chain A, published generic humanized VHH scaffold','framework_only_ids':[r['id'] for r in framework],'deleted_regions_in_source_numbering':regions,'loop_initialization':'independent random non-antibody backbone torsion fragments; no donor CDR sequence or coordinates','framework_FERA_option':'47=E,48=R; hydrophilic VHH hallmark alternative','terminal_append':'S'},indent=2))
fxyz=np.array([p for r in framework for name,p in r['atoms'].items() if name!='OXT']);fids=np.array([r['id'] for r in framework for name in r['atoms'] if name!='OXT'])
frags=[]
for name in ['human6ARU','1IVO','1NQL','1YY9','2A3D','1FSD','1L2Y']:
 rr=json.loads((O/(name+'_chainA.json')).read_text());vals=[]
 for j in range(1,len(rr)-1):
  try:
   prev=rr[j-1]['atoms'];cur=rr[j]['atoms'];nex=rr[j+1]['atoms']
   if np.linalg.norm(np.array(prev['C'])-cur['N'])>1.6 or np.linalg.norm(np.array(cur['C'])-nex['N'])>1.6:vals.append(None);continue
   phi=dihedral(np.array(prev['C']),np.array(cur['N']),np.array(cur['CA']),np.array(cur['C']));psi=dihedral(np.array(cur['N']),np.array(cur['CA']),np.array(cur['C']),np.array(nex['N']))
   vals.append([phi,psi])
  except KeyError:vals.append(None)
 for i in range(len(vals)-2):
  if all(v is not None for v in vals[i:i+3]):frags.append(vals[i:i+3])
frags=np.array(frags);np.save(R/'stage2/intermediate/nonbinder_torsion_fragments.npy',frags)
a=np.array([0.,0.,0.]);b=np.array([1.5,0.,0.]);c=np.array([2.,1.,0.])
for phi in [-2.,-.5,1.,2.5]:
 d=place(a,b,c,1.5,1.9,phi)
 assert abs(((dihedral(a,b,c,d)-phi+np.pi)%(2*np.pi)-np.pi))<1e-8
rng=np.random.Generator(np.random.MT19937(20031998))
configs=[('H1',11),('H1',13),('H2',13),('H2',15),('H3',14),('H3',16),('H3',18),('H3',20)]
tries=int(sys.argv[1]) if len(sys.argv)>1 else 1400
wanted=int(sys.argv[2]) if len(sys.argv)>2 else 35
stats={}
for reg,n in configs:
 lo,hi=regions[reg];l=byid[lo-1]['atoms'];r=byid[hi+1]['atoms'];anchor=np.array([l[k] for k in ['N','CA','C']]);target=np.array([r[k] for k in ['N','CA','C']])
 psi_anchor=dihedral(anchor[0],anchor[1],anchor[2],np.array(l['O']))-np.pi
 keep=[];counts={'attempts':0,'closed':0,'rama_ok':0,'steric_ok':0};t=time.time()
 for k in range(tries):
  pp=np.concatenate([frags[rng.integers(len(frags))] for _ in range((n+3)//3)],axis=0)[:n+1]
  pp=pp+rng.normal(0,.12,pp.shape)
  x=make_chain(anchor,pp[:,0],pp[:-1,1],psi_anchor,target)
  x,err,cycles=ccd(x,target)
  counts['attempts']+=1
  if err>.09:continue
  counts['closed']+=1;bb,angles=backbone_atoms(x);q,pos=rama_quality(angles)
  if q.max()>12 or q.mean()>2.6 or pos>max(4,n*.35):continue
  counts['rama_ok']+=1
  bad,mind=loop_clashes(bb,fxyz,fids,lo-1,hi+1)
  if bad:continue
  counts['steric_ok']+=1
  if keep and min(np.sqrt(np.mean(np.sum((bb[:,1]-old['bb'][:,1])**2,axis=1))) for old in keep)<.65:continue
  keep.append({'bb':bb,'pp':angles,'rama':q,'closure':err,'cycles':cycles,'source_attempt':k})
  if len(keep)>=wanted:break
 print(reg,n,counts,'accepted',len(keep),'seconds',round(time.time()-t,2),flush=True)
 stats[reg+str(n)]={**counts,'accepted':len(keep),'seconds':time.time()-t}
 if keep:
  np.savez_compressed(out/(reg+str(n)+'.npz'),bb=np.array([a['bb'] for a in keep]),pp=np.array([a['pp'] for a in keep]),rama=np.array([a['rama'] for a in keep]),closure=np.array([a['closure'] for a in keep]),attempt=np.array([a['source_attempt'] for a in keep]))
(R/'stage2/reference/loop_generation_stats.json').write_text(json.dumps(stats,indent=2))
