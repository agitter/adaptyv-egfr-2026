"""Framework-aligned end-point contact audit of unbound thermal stress tests.
Least-squares rigid alignment (Kabsch1976); no relaxation or affinity inference.
"""
from pathlib import Path
import sys,json,hashlib
import numpy as np
from scipy.spatial import cKDTree
R=Path('/mnt/data/egfr_campaign');S=R/'stage4'
sys.path[:0]=[str(R/'stage2/code'),str(R/'stage2/runtime/python')]
from mm_refine import app,u

def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def atoms(path,chain):
 p=app.PDBFile(str(path));x=np.asarray(p.positions.value_in_unit(u.angstrom));return [{a.name:x[a.index] for a in r.atoms() if a.element!=app.element.hydrogen} for r in p.topology.residues() if chain is None or r.chain.id==chain]
def align(moving,fixed):
 xc=moving.mean(0);yc=fixed.mean(0);U,ss,V=np.linalg.svd((moving-xc).T@(fixed-yc));z=np.eye(3);z[-1,-1]=np.linalg.det(U@V);rot=U@z@V;return rot,yc-xc@rot

def run(tag):
 traj=S/'intermediate/dynamics'/(tag+'.json');d=json.loads(traj.read_text());cid=d['candidate_id']
 if cid=='NATIVE_3EAK':return
 root=next(r for r in [S,R/'stage3'] if (r/'intermediate/designs'/(cid+'.json')).exists());design=json.loads((root/'intermediate/designs'/(cid+'.json')).read_text());source=Path(d['metadata']['source']);final=S/'intermediate/dynamics'/(tag+'_final.pdb');ref=atoms(source,'B');free=atoms(final,None);target=atoms(source,'A');assert len(ref)==len(free)==len(design['sequence'])
 framework=np.array([not r['loop'] for r in design['structure']]);rot,trans=align(np.array([a['CA'] for a in free])[framework],np.array([a['CA'] for a in ref])[framework]);aligned=[{n:x@rot+trans for n,x in r.items()} for r in free]
 human=json.loads((R/'intermediate/human6ARU_aligned.json').read_text());alltarget=[]
 for r in human:
  ar=target[r['human_pos']-334] if 334<=r['human_pos']<=505 else r['atoms']
  alltarget.extend(list(ar.values()))
 gly=json.loads((R/'stage2/intermediate/target_glycans.json').read_text());glyxyz=np.array([g['xyz'] for g in gly]);targetxyz=np.array(alltarget);rows=[]
 for pos in [52,54,103]:
  aa=design['sequence'][pos-1]
  if aa not in 'DE':continue
  names=['OD1','OD2'] if aa=='D' else ['OE1','OE2']
  for hp in range(334,506):
   if not 334<=hp<=505:continue
   rr=target[hp-334]
   if not {'ND1','NE2'}<=set(rr):continue
   h=np.array([rr[n] for n in ['ND1','NE2']]);old=np.array([ref[pos-1][n] for n in names]);new=np.array([aligned[pos-1][n] for n in names]);a=float(np.linalg.norm(old[:,None,:]-h[None,:,:],axis=2).min());b=float(np.linalg.norm(new[:,None,:]-h[None,:,:],axis=2).min())
   rows.append({'binder_position':pos,'binder_residue':aa,'target_human_position':hp,'initial_min_ON_A':a,'final_min_ON_A':b,'initial_contact_below4A':a<4.,'final_contact_below4A':b<4.})
 finalxyz=np.array([x for rr in aligned for x in rr.values()]);distance=cKDTree(targetxyz).query(finalxyz)[0];gd=cKDTree(glyxyz).query(finalxyz)[0]
 result={'trajectory':tag,'candidate_id':cid,'source_sha256':digest(source),'final_sha256':digest(final),'framework_aligned':True,'contacts':rows,'initial_contacts_below4A':sum(r['initial_contact_below4A'] for r in rows),'retained_initial_contacts_below4A':sum(r['initial_contact_below4A'] and r['final_contact_below4A'] for r in rows),'final_receptor_heavy_atoms_below2A':int(sum(distance<2.)),'final_min_receptor_A':float(distance.min()),'final_glycan_heavy_atoms_below2A':int(sum(gd<2.)),'final_min_glycan_A':float(gd.min()),'qualification':'End-point geometry after a very short unbound trajectory; no target relaxation. Not a binding, kinetic or equilibrium-stability prediction.'}
 out=S/'intermediate/trajectory_contacts';out.mkdir(exist_ok=True);(out/(tag+'.json')).write_text(json.dumps(result,indent=2));print(tag,result['retained_initial_contacts_below4A'],'/',result['initial_contacts_below4A'],'target clashes',result['final_receptor_heavy_atoms_below2A'],flush=True);return result

def tests():
 rng=np.random.Generator(np.random.MT19937(1976));a=rng.normal(size=(20,3));theta=.4;rot=np.array([[np.cos(theta),-np.sin(theta),0],[np.sin(theta),np.cos(theta),0],[0,0,1]]);b=a@rot+np.array([4,-7,2]);rr,tt=align(a,b);err=float(np.max(abs(a@rr+tt-b)));assert err<1e-10;return {'proper_rotation_alignment_max_error_A':err,'pass':True}
if __name__=='__main__':
 if len(sys.argv)>1:run(sys.argv[1])
 else:
  (S/'reference/trajectory_contact_tests.json').write_text(json.dumps(tests(),indent=2))
  for p in (S/'intermediate/dynamics').glob('*_v2_s*.json'):
   if not p.name.endswith('_progress.json'):run(p.stem)
