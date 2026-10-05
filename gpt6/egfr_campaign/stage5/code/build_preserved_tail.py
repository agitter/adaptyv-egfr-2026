"""Own-design loop returns: CCD2003 + non-antibody torsion fragments.
No native antibody loop coordinates/sequences initialize the rebuilt atoms.
"""
from pathlib import Path
import sys,json,math,copy,hashlib,time
import numpy as np
from scipy.spatial import cKDTree
R=Path('/mnt/data/egfr_campaign');S5=R/'stage5'
sys.path[:0]=[str(R/'stage2/code'),str(R/'code'),str(R/'stage2/runtime/python')]
from legacy_geometry import make_chain,ccd,backbone_atoms,rama_quality,loop_clashes,dihedral
from mm_refine import app,u
PARENT=R/'stage3/intermediate/designs/C00003.json'
SOURCE=R/'stage4/intermediate/matched_controls/intermediate/refined/C00003_human6ARU_relaxed.pdb'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def load_parent():
 d=json.loads(PARENT.read_text());p=app.PDBFile(str(SOURCE));x=np.asarray(p.positions.value_in_unit(u.angstrom));rr=[r for r in p.topology.residues() if r.chain.id=='B'];rows=copy.deepcopy(d['structure'])
 for row,r in zip(rows,rr):
  row['atoms']={a.name:x[a.index].tolist() for a in r.atoms() if a.element!=app.element.hydrogen and a.name!='OXT'};row['parent_position']=row['position']
 target=[x[a.index] for r in p.topology.residues() if r.chain.id=='A' for a in r.atoms() if a.element!=app.element.hydrogen]
 return d,rows,np.array(target)
def run(region,n,tries=1200,wanted=12):
 t=time.time();out=S5/'intermediate/backbones';out.mkdir(parents=True,exist_ok=True);tag=f'{region}_n{n:02d}';dest=out/(tag+'.npz')
 if dest.exists():raise FileExistsError(dest)
 d,rows,tx=load_parent();left,right={'T3':(107,116)}[region]
 l=rows[left-1]['atoms'];r=rows[right-1]['atoms'];anchor=np.array([l[k] for k in ['N','CA','C']]);target=np.array([r[k] for k in ['N','CA','C']]);psi=dihedral(*anchor,np.array(l['O']))-math.pi
 fxyz=[];fids=[]
 for idx,row in enumerate(rows,1):
  if left<idx<right:continue
  for name,p in row['atoms'].items():fxyz.append(p);fids.append(idx)
 fxyz=np.array(fxyz);fids=np.array(fids)
 mouse=json.loads((R/'intermediate/mouseAF_aligned.json').read_text());mx=np.array([p for row in mouse for p in row['atoms'].values()]);gly=json.loads((R/'stage2/intermediate/target_glycans.json').read_text());gx=np.array([a['xyz'] for a in gly]);trees=[cKDTree(x) for x in [tx,mx,gx]]
 fragments=np.load(R/'stage2/intermediate/nonbinder_torsion_fragments.npy');seed=200319980+100*4+n;rng=np.random.Generator(np.random.MT19937(seed));kept=[];stats={'attempts':0,'closed':0,'rama_pass':0,'self_steric_pass':0,'target_glycan_pass':0,'diverse':0};meta=[]
 for trial in range(tries):
  pp=np.concatenate([fragments[int(rng.integers(len(fragments)))] for _ in range((n+3)//3)])[:n+1]+rng.normal(0,.08,(n+1,2))
  x=make_chain(anchor,pp[:,0],pp[:-1,1],psi,target);x,err,cycles=ccd(x,target,600,.065);stats['attempts']+=1
  if err>.075:continue
  stats['closed']+=1;bb,angles=backbone_atoms(x);q,pos=rama_quality(angles)
  if q.max()>8.0 or q.mean()>2.0 or pos>max(2,n*.28):continue
  stats['rama_pass']+=1;clashes,mind=loop_clashes(bb,fxyz,fids,left,right)
  if clashes:continue
  stats['self_steric_pass']+=1;bx=bb[:,:4,:].reshape(-1,3);mins=[float(tr.query(bx)[0].min()) for tr in trees]
  if min(mins[:2])<2.15 or mins[2]<2.6:continue
  stats['target_glycan_pass']+=1
  if kept and min(float(np.sqrt(np.mean(np.sum((bb[:,1]-b[:,1])**2,1)))) for b in kept)<.75:continue
  kept.append(bb);meta.append({'source_attempt':trial,'closure_A':float(err),'cycles':int(cycles),'rama':q.tolist(),'angles_radians':angles.tolist(),'positive_phi_count':int(pos),'min_backbone_distances_A':dict(zip(['human','mouse','resolved_glycans'],mins))});stats['diverse']+=1
  if len(kept)>=wanted:break
 if kept:np.savez_compressed(dest,bb=np.array(kept))
 info={'region':region,'left_parent_anchor':left,'right_parent_anchor':right,'new_segment_length':n,'parent_segment_length':right-left-1,'tag':tag,'seed':seed,'rng':'MT19937','stats':stats,'accepted':meta,'seconds':time.time()-t,'parent_sequence_file':str(PARENT),'parent_sequence_sha256':sha(PARENT),'parent_geometry_file':str(SOURCE),'parent_geometry_sha256':sha(SOURCE),'torsion_fragment_file':str(R/'stage2/intermediate/nonbinder_torsion_fragments.npy'),'torsion_fragment_sha256':sha(R/'stage2/intermediate/nonbinder_torsion_fragments.npy'),'method':'CCD2003; non-antibody phi/psi fragments; own-design fixed anchors; no native antibody CDR'}
 (out/(tag+'.json')).write_text(json.dumps(info,indent=2));print(tag,stats,'seconds',round(info['seconds'],1),flush=True);return info
if __name__=='__main__':run(sys.argv[1],int(sys.argv[2]),int(sys.argv[3]) if len(sys.argv)>3 else 1000,int(sys.argv[4]) if len(sys.argv)>4 else 18)
