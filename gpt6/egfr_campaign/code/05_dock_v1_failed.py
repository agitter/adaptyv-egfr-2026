"""Classical random rigid-body docking of rebuilt generic VHH frameworks.
Backbone library and rotations are sampled independently of known binders.
"""
from pathlib import Path
import json,numpy as np,time,sys
from scipy.spatial import cKDTree
from scipy.spatial.transform import Rotation
from legacy_geometry import unit,atom_frame
R=Path('/mnt/data/egfr_campaign');O=R/'intermediate';OUT=O/'docking';OUT.mkdir(exist_ok=True)
rng=np.random.Generator(np.random.MT19937(19762003))
framework=json.loads((O/'framework_only.json').read_text());fd={x['id']:x for x in framework}
regions={'H1':(26,38),'H2':(53,67),'H3':(100,117)}
libs={}
for p in (O/'loops').glob('*.npz'):
 with np.load(p) as z:libs[p.stem]={k:z[k] for k in z.files}
choices={reg:[k for k in libs if k.startswith(reg)] for reg in regions}
targets={name:json.loads((O/(name+'_aligned.json')).read_text()) for name in ['human6ARU','mouseAF','1NQL','1IVO']}
ht=targets['human6ARU'];hd={x['human_pos']:x for x in ht}
xyz={n:np.array([v for r in a for k,v in r['atoms'].items()]) for n,a in targets.items()};trees={n:cKDTree(v) for n,v in xyz.items()}
thpos=np.array([r['human_pos'] for r in ht for k in r['atoms']]);thcons=np.array([r['conserved'] for r in ht for k in r['atoms']])
gly=json.loads((O/'human6ARU_glycans.json').read_text());glyxyz=np.array([g['xyz'] for g in gly]);gtree=cKDTree(glyxyz)
# Protect potential glycans at sequons, not only crystallographically resolved sugars.
sequons=json.loads((R/'reference/glycosylation_sequons.json').read_text())['human'];gn=[]
for pos in sequons:
 if pos in hd and 'ND2' in hd[pos]['atoms']:gn.append(hd[pos]['atoms']['ND2'])
gn=np.array(gn);gntree=cKDTree(gn)
patches=[('A',[370,433,368]),('B',[433,460,458]),('C',[370,433]),('D',[358,344,347]),('E',[433,489,460]),('F',[344,347,370])]
charges={}
for r in ht:
 if r['aa'] in 'HDE' and 334<=r['human_pos']<=505 and r.get('sasa',0)>15:
  names={'H':['ND1','NE2'],'D':['OD1','OD2'],'E':['OE1','OE2']}[r['aa']]
  if all(k in r['atoms'] for k in names):charges[r['human_pos']]=np.mean([r['atoms'][k] for k in names],axis=0)
cc=np.array([charges[p] for p in charges]);cpos=list(charges)
patchinfo=[]
for name,positions in patches:
 point=np.mean([charges[p] if p in charges else hd[p]['atoms']['CA'] for p in positions],axis=0)
 near=np.array([r['atoms']['CA'] for r in ht if 334<=r['human_pos']<=505 and np.linalg.norm(np.array(r['atoms']['CA'])-point)<17])
 z=unit(point-near.mean(0));x=unit(np.cross(z,np.array([0.,0.,1.])))
 if np.linalg.norm(x)<.5:x=unit(np.cross(z,np.array([0.,1.,0.])))
 y=np.cross(z,x);patchinfo.append((name,point,np.array([x,y,z])))

def assemble(sel):
 bb=[];mask=[];aa=[];orig=[];phipsi=[];frames=[];qsum=0.;clerr=0.
 for lo,hi,key in [(1,25,None),(26,38,'H1'),(39,52,None),(53,67,'H2'),(68,99,None),(100,117,'H3'),(118,127,None)]:
  if key:
   lib,idx=sel[key];b=libs[lib]['bb'][idx];pp=libs[lib]['pp'][idx]
   for i in range(len(b)):
    bb.append(b[i]);mask.append(True);aa.append('X');orig.append(1000+len(bb));phipsi.append(pp[i])
   qsum+=float(libs[lib]['rama'][idx].sum());clerr=max(clerr,float(libs[lib]['closure'][idx]))
  else:
   for i in range(lo,hi+1):
    a=fd[i]['atoms'];row=[a[k] for k in ['N','CA','C','O']]
    row.append(a.get('CB',a['CA']));bb.append(row);mask.append(False);aa.append(fd[i]['aa']);orig.append(i);phipsi.append([0.,0.])
    for n,v in a.items():frames.append((len(bb)-1,n,v))
 return np.array(bb),np.array(mask),''.join(aa),np.array(orig),np.array(phipsi),frames,qsum/sum(mask),clerr

def dock_score(bb,mask):
 points=bb[:,:4].reshape(-1,3)
 dh,idx=trees['human6ARU'].query(points,k=1)
 if dh.min()<2.5:return None
 # No design can fit human by intersecting the mouse ECR.
 dm,_=trees['mouseAF'].query(points,k=1)
 if dm.min()<2.25:return None
 dg,_=gtree.query(points,k=1)
 if dg.min()<2.5:return None
 glydist,_=gntree.query(bb[mask,1],k=1)
 if glydist.min()<4.5:return None
 cp=bb[mask,4];ca=bb[mask,1];dire=cp-ca; dire/=np.maximum(np.linalg.norm(dire,axis=1)[:,None],1e-8)
 d,ix=trees['human6ARU'].query(cp,k=1)
 ncontact=np.sum(d<5.8)
 if ncontact<6 or ncontact>28:return None
 diff=cc[None,:,:]-cp[:,None,:];dd=np.linalg.norm(diff,axis=2)
 orientation=np.sum(diff*dire[:,None,:],axis=2)/np.maximum(dd,1e-6)
 opp=np.exp(-((dd-4.6)/2.8)**2)*np.clip(orientation+.4,0,1)
 best=opp.max(0)
 hscore=sum(best[k] for k,p in enumerate(cpos) if hd[p]['aa']=='H' and hd[p]['conserved'])
 ascore=sum(sorted([best[k] for k,p in enumerate(cpos) if hd[p]['aa'] in 'DE' and hd[p]['conserved']],reverse=True)[:4])
 # Reward independent titratable opportunities, rather than histidine abundance.
 if hscore<.30 and ascore<1.4:return None
 gap=np.sum(np.clip(3.-dh,0,None)**2)*15+np.sum(np.clip(2.8-dm,0,None)**2)*10
 noncons=np.sum((d<5.3)*(~thcons[ix]));framework_contact=np.sum(dh.reshape(-1,4)[~mask].min(1)<4.)
 score=3.0*hscore+1.2*ascore+.13*ncontact-.9*noncons-.35*framework_contact-gap-.2*np.sum(glydist<7.)
 return dict(score=float(score),h_opportunity=float(hscore),acid_opportunity=float(ascore),contact_residues=int(ncontact),nonconserved_nearest=int(noncons),minimum_human_backbone_distance=float(dh.min()),minimum_mouse_backbone_distance=float(dm.min()),glycan_anchor_min=float(glydist.min()))

N=int(sys.argv[1]) if len(sys.argv)>1 else 500
TRIALS=int(sys.argv[2]) if len(sys.argv)>2 else 240
results=[];attempts=0;accepted=0;t0=time.time()
while accepted<N and attempts<N*15:
 attempts+=1;sel={reg:(rng.choice(choices[reg]),None) for reg in regions}
 sel={r:(str(k),int(rng.integers(len(libs[k]['bb'])))) for r,(k,_) in sel.items()}
 bb,mask,aas,orig,pp,frames,rama,closure=assemble(sel)
 # Check all rebuilt loops against one another before any target docking.
 loops=[np.where(mask)[0]];n=len(bb);i,j=np.triu_indices(n,k=2)
 where=mask[i]&mask[j]
 ii=i[where];jj=j[where]
 if np.any(np.linalg.norm(bb[ii,1]-bb[jj,1],axis=1)<3.2):continue
 # Atom-pair checks of nonlocal loop backbones.
 bad=False
 for reg in regions:
  pass
 mut=np.flatnonzero(mask)
 for a in range(len(mut)):
  i=mut[a]
  others=mut[a+1:];others=others[others-i>1]
  if not len(others):continue
  dd=np.linalg.norm(bb[i,:4,None,None,:]-bb[others,None,:4,:],axis=-1)
  if dd.min()<2.15:bad=True;break
 if bad:continue
 accepted+=1
 center=bb[mask,1].mean(0);fc=bb[~mask,1].mean(0);z=unit(center-fc);x=unit(bb[25,1]-center-z*np.dot(bb[25,1]-center,z));y=np.cross(z,x);basis=np.array([x,y,z])
 best=[]
 for trial in range(TRIALS):
  patch,point,tbasis=patchinfo[trial%len(patchinfo)]
  # Orient the paratope toward the receptor and scan azimuth/tilt/offset.
  twist=rng.uniform(-np.pi,np.pi);tilt=rng.normal(0,.23,2)
  rr=Rotation.from_rotvec(np.array([tilt[0],tilt[1],twist])).as_matrix()
  flip=np.diag([1.,-1.,-1.]);rot=basis.T@rr.T@flip@tbasis
  tr=point+tbasis[2]*rng.uniform(6.,10.)+tbasis[0]*rng.normal(0,2.7)+tbasis[1]*rng.normal(0,2.7)-center@rot
  b=bb@rot+tr
  metrics=dock_score(b,mask)
  if metrics is None:continue
  metrics['score']-=.15*rama
  rec={'backbone':accepted,'patch':patch,'selection':sel,'rotation':rot.tolist(),'translation':tr.tolist(),'rama_mean':rama,'closure_max':closure,**metrics}
  best.append(rec)
 if best:
  best.sort(key=lambda x:x['score'],reverse=True)
  # Keep up to three poses with different geometry from each backbone.
  chosen=[]
  for rec in best:
   r=np.array(rec['rotation']);tr=np.array(rec['translation']);pts=bb[mask,1]@r+tr
   if all(np.sqrt(np.mean((pts-prev)**2))>2. for prev in chosen):
    chosen.append(pts);results.append(rec)
   if len(chosen)==3:break
 if accepted%25==0:
  print('assembled',accepted,'attempts',attempts,'poses',len(results),'elapsed',round(time.time()-t0,1),flush=True)
  (OUT/'poses_checkpoint.json').write_text(json.dumps(sorted(results,key=lambda x:x['score'],reverse=True)))
results.sort(key=lambda x:x['score'],reverse=True)
for i,r in enumerate(results):r['pose_id']='P%05d'%(i+1)
(OUT/'poses.json').write_text(json.dumps(results))
(R/'reference/docking_stats.json').write_text(json.dumps({'backbones_assembled':accepted,'assembly_attempts':attempts,'orientations_per_backbone':TRIALS,'poses_retained':len(results),'seconds':time.time()-t0,'seed':19762003,'rng':'MT19937'},indent=2))
print('DONE',accepted,len(results),round(time.time()-t0,1),flush=True)
