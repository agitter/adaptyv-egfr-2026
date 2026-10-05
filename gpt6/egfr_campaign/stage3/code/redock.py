"""Rigid-body receptor-His-focused docking, using the stage-2 MT19937 loops.
All original antigen-binding loops were removed before those loops were made.
No known binder orientation or binding-loop coordinates enter this search.
"""
from pathlib import Path
import sys,json,time,math
import numpy as np
from scipy.spatial import cKDTree
from scipy.spatial.transform import Rotation
R=Path('/mnt/data/egfr_campaign');S=R/'stage3'
sys.path[:0]=[str(R/'stage2/code'),str(R/'code')]
import dock_mt as d
from classical_design import charge_center
from legacy_geometry import unit


def main():
    rng=np.random.Generator(np.random.MT19937(196419982003))
    records=json.loads((R/'stage2/intermediate/docking/poses.json').read_text())
    gly=json.loads((R/'stage2/intermediate/target_glycans.json').read_text());gt=cKDTree([x['xyz'] for x in gly])
    outside=cKDTree([p for r in d.targets['1NQL'] if not 334<=r['human_pos']<=505 for p in r['atoms'].values()])
    hpos=[358,370,383,433];hc=np.array([charge_center('H',d.hd[p]['atoms']) for p in hpos])
    patches=[]
    for label,sites in [('G',[358,370,383]),('H',[370,383,433]),('I',[358,383]),('J',[370,433]),('K',[358,370])]:
        point=np.mean([charge_center('H',d.hd[p]['atoms']) for p in sites],axis=0)
        near=np.array([r['atoms']['CA'] for r in d.ht if 334<=r['human_pos']<=505 and np.linalg.norm(np.array(r['atoms']['CA'])-point)<18])
        z=unit(point-near.mean(0));x=unit(np.cross(z,[0.,0.,1.]));y=np.cross(z,x);patches.append((label,point,np.array([x,y,z])))
    seen=set();unique=[]
    for rec in records:
        key=json.dumps(rec['selection'],sort_keys=True)
        if key not in seen:unique.append(rec);seen.add(key)
    # Retain diverse original backbones; their rank does not dictate the new orientation.
    order=rng.permutation(len(unique))[:240];output=[];attempts=0;passed=0;start=time.time()
    for ix in order:
        src=unique[ix];bb,mask,aa,orig,pp,frames,rama,closure=d.assemble(src['selection'])
        fixed=np.array([p for _,_,p in frames]);center=bb[mask,1].mean(0);fc=bb[~mask,1].mean(0)
        z=unit(center-fc);x=unit(bb[25,1]-center-z*np.dot(bb[25,1]-center,z));basis=np.array([x,np.cross(z,x),z])
        local=[]
        for trial in range(360):
            attempts+=1;label,point,tbasis=patches[trial%len(patches)]
            twist=rng.uniform(-np.pi,np.pi);tilt=rng.normal(0,.28,2)
            rr=Rotation.from_rotvec([tilt[0],tilt[1],twist]).as_matrix()
            rot=basis.T@rr.T@np.diag([1.,-1.,-1.])@tbasis
            tr=point+tbasis[2]*rng.uniform(5.,10.)+tbasis[0]*rng.normal(0,3.2)+tbasis[1]*rng.normal(0,3.2)-center@rot
            coords=bb@rot+tr;metrics=d.dock_score(coords,mask)
            if metrics is None:continue
            allpoints=np.concatenate([coords[:,:4].reshape(-1,3),fixed@rot+tr])
            if gt.query(allpoints)[0].min()<2.45 or outside.query(allpoints)[0].min()<2.1:continue
            if min(d.trees[n].query(fixed@rot+tr)[0].min() for n in ['human6ARU','mouseAF'])<2.2:continue
            cb=coords[mask,4];ca=coords[mask,1];out=(cb-ca)/np.maximum(np.linalg.norm(cb-ca,axis=1)[:,None],1e-9)
            delta=hc[None,:,:]-cb[:,None,:];dist=np.linalg.norm(delta,axis=2)
            orient=np.sum(delta*out[:,None,:],axis=2)/np.maximum(dist,1e-9)
            opportunity=(np.exp(-((dist-4.8)/2.8)**2)*np.clip(orient+.3,0,1)).max(0)
            nh=int((opportunity>.4).sum())
            if nh<2:continue
            passed+=1
            score=float(3.0*opportunity.sum()+.10*metrics['contact_residues']-.4*metrics['nonconserved_nearest']-.15*rama)
            rec={'selection':src['selection'],'rotation':rot.tolist(),'translation':tr.tolist(),'patch':label,
                 'rama_mean':rama,'closure_max':closure,**metrics,'score':score,'histidine_positions':hpos,
                 'histidine_opportunity':opportunity.tolist(),'independent_his_opportunities':nh,
                 'ancestor_backbone_pose':src['pose_id'],'generation':'stage3 receptor-His rigid docking'}
            local.append((rec,coords[mask,1]))
        local.sort(key=lambda v:v[0]['score'],reverse=True);chosen=[]
        for rec,xyz in local:
            if all(np.sqrt(np.mean(np.sum((xyz-other)**2,axis=1)))>3.0 for other in chosen):output.append(rec);chosen.append(xyz)
            if len(chosen)>=3:break
        if len(output)%50<3:print('orientations',attempts,'retained',len(output),'seconds',round(time.time()-start,1),flush=True)
    # Original poses remain a separate, traceable branch and undergo the new filters.
    for rec in records:
        r=rec.copy();r['ancestor_backbone_pose']=rec['pose_id'];r['generation']='stage2 MT19937 pose, redesigned from empty CDR sequence';output.append(r)
    output.sort(key=lambda r:r['score'],reverse=True)
    for n,rec in enumerate(output,1):rec['source_pose_id']=rec.get('pose_id');rec['pose_id']='G%05d'%n
    dest=S/'intermediate/docking';dest.mkdir(parents=True,exist_ok=True)
    (dest/'poses.json').write_text(json.dumps(output))
    stats={'orientations':attempts,'passing_orientations':passed,'new_poses':len(output)-len(records),'inherited_poses':len(records),'seed':196419982003,'rng':'MT19937','seconds':time.time()-start}
    (S/'reference/redocking.json').write_text(json.dumps(stats,indent=2));print(stats,flush=True)

if __name__=='__main__':main()
