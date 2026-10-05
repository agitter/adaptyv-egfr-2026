"""Classical ideal-beta return extension plus CCD loop closure.
No native antibody CDR atom initializes any new residue. Standard peptide
geometry/Ramachandran beta region, reverse internal coordinates and CCD2003.
This is not a claim of a canonical antibody loop or a folded beta strand.
"""
from pathlib import Path
import sys,math,json,time
import numpy as np
from scipy.spatial import cKDTree
S=Path('/mnt/data/egfr_campaign/stage5');R=S.parent;sys.path[:0]=[str(S/'code'),str(R/'stage2/code'),str(R/'code'),str(R/'stage2/runtime/python')]
from build_return_loops import load_parent,PARENT,SOURCE,sha
from legacy_geometry import place,dihedral,make_chain,ccd,backbone_atoms,rama_quality,loop_clashes

def extend_backwards(anchor,phi_right,phi_beta,psi_beta,count=3):
 current={k:np.array(v) for k,v in anchor.items()};parts=[]
 for j in range(count):
  ph=phi_right if j==0 else phi_beta
  C=place(current['C'],current['CA'],current['N'],1.329,math.radians(121.7),ph)
  CA=place(current['CA'],current['N'],C,1.525,math.radians(116.2),math.pi)
  N=place(current['N'],C,CA,1.458,math.radians(111.2),psi_beta)
  O=place(N,CA,C,1.231,math.radians(120.8),psi_beta+math.pi)
  row={'N':N,'CA':CA,'C':C,'O':O};parts.insert(0,row);current=row
 return parts

def run(n,tries=1600,wanted=16):
 start=time.time();d,rows,tx=load_parent();left,right=54,66;leftbb=np.array([rows[left-1]['atoms'][k] for k in ['N','CA','C']]);rightrow=rows[right-1]['atoms'];rightphi=dihedral(np.array(rows[right-2]['atoms']['C']),*[np.array(rightrow[k]) for k in ['N','CA','C']]);psi=dihedral(*leftbb,np.array(rows[left-1]['atoms']['O']))-math.pi
 fx=[];fi=[]
 for i,row in enumerate(rows,1):
  if left<i<right:continue
  for v in row['atoms'].values():fx.append(v);fi.append(i)
 fx=np.array(fx);fi=np.array(fi);mouse=json.loads((R/'intermediate/mouseAF_aligned.json').read_text());mx=np.array([v for r in mouse for v in r['atoms'].values()]);gx=np.array([g['xyz'] for g in json.loads((R/'stage2/intermediate/target_glycans.json').read_text())]);trees=[cKDTree(v) for v in [tx,mx,gx]]
 fragments=np.load(R/'stage2/intermediate/nonbinder_torsion_fragments.npy');seed=1983200300+n;rng=np.random.Generator(np.random.MT19937(seed));kept=[];meta=[];stats={'attempts':0,'closed':0,'rama_pass':0,'self_steric_pass':0,'target_glycan_pass':0,'diverse':0};nflex=n-3
 for trial in range(tries):
  # Standard beta region, not coordinates sampled from native antibody loops.
  ph=math.radians(float(rng.uniform(-145,-110)));ps=math.radians(float(rng.uniform(115,150)));ext=extend_backwards(rightrow,rightphi,ph,ps,3);target=np.array([ext[0][k] for k in ['N','CA','C']]);pp=np.concatenate([fragments[int(rng.integers(len(fragments)))] for _ in range((nflex+3)//3)])[:nflex+1]+rng.normal(0,.08,(nflex+1,2));x=make_chain(leftbb,pp[:,0],pp[:-1,1],psi,target);x,err,cycles=ccd(x,target,600,.065);stats['attempts']+=1
  if err>.075:continue
  stats['closed']+=1;flexbb,angles=backbone_atoms(x);fullbb=np.zeros((n,5,3));fullbb[:nflex]=flexbb
  for j,row in enumerate(ext):
   for k,name in enumerate(['N','CA','C','O']):fullbb[nflex+j,k]=row[name]
  allangles=[]
  for i in range(n):
   prev=leftbb[2] if i==0 else fullbb[i-1,2];nxt=np.array(rightrow['N']) if i==n-1 else fullbb[i+1,0];allangles.append([dihedral(prev,*fullbb[i,:3]),dihedral(*fullbb[i,:3],nxt)])
  allangles=np.array(allangles);q,pos=rama_quality(allangles)
  if q.max()>8 or q.mean()>2 or pos>max(2,n*.28):continue
  stats['rama_pass']+=1;bad,mind=loop_clashes(fullbb,fx,fi,left,right)
  if bad:continue
  stats['self_steric_pass']+=1;bx=fullbb[:,:4].reshape(-1,3);mins=[float(t.query(bx)[0].min()) for t in trees]
  if min(mins[:2])<2.15 or mins[2]<2.6:continue
  stats['target_glycan_pass']+=1
  if kept and min(float(np.sqrt(np.mean(np.sum((fullbb[:,1]-b[:,1])**2,1)))) for b in kept)<.65:continue
  kept.append(fullbb);meta.append({'source_attempt':trial,'closure_A':float(err),'cycles':int(cycles),'rama':q.tolist(),'angles_radians':allangles.tolist(),'positive_phi_count':int(pos),'ideal_extension_phi_deg':math.degrees(ph),'ideal_extension_psi_deg':math.degrees(ps),'retained_right_framework_phi_deg':math.degrees(rightphi),'min_backbone_distances_A':dict(zip(['human','mouse','resolved_glycans'],mins))});stats['diverse']+=1
  if len(kept)>=wanted:break
 out=S/'intermediate/backbones';tag=f'B2_n{n:02d}';dest=out/(tag+'.npz')
 if dest.exists():raise FileExistsError(dest)
 if kept:np.savez_compressed(dest,bb=np.array(kept))
 info={'region':'H2','left_parent_anchor':left,'right_parent_anchor':right,'new_segment_length':n,'parent_segment_length':11,'tag':tag,'seed':seed,'rng':'MT19937','stats':stats,'accepted':meta,'seconds':time.time()-start,'parent_sequence_file':str(PARENT),'parent_sequence_sha256':sha(PARENT),'parent_geometry_file':str(SOURCE),'parent_geometry_sha256':sha(SOURCE),'method':'Standard ideal beta peptide geometry plus independently generated non-antibody torsion CCD2003; no native CDR copying; fold unvalidated.'};(out/(tag+'.json')).write_text(json.dumps(info,indent=2));print(tag,stats,round(info['seconds'],1),flush=True);return info
if __name__=='__main__':run(int(sys.argv[1]))
