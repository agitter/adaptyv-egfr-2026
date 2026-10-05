"""Classical discrete rotamer packing of newly closed loop-return segments.
Metropolis/annealing (1953/1983), empirical rotamers, distance/H-bond/solvent
heuristics. Coarse scores are filters, not free energies. No sequence models.
"""
from pathlib import Path
import sys,json,copy,math,re,time,hashlib
import numpy as np
from numba import njit
from scipy.spatial import cKDTree
R=Path('/mnt/data/egfr_campaign');S5=R/'stage5'
sys.path[:0]=[str(S5/'code'),str(R/'stage2/code'),str(R/'code'),str(R/'stage2/runtime/python')]
from build_return_loops import load_parent,PARENT,SOURCE,sha
from classical_design import LIB,BB,conformer,atom_data,pair_score,charge_center,screened_charge,CHARGES
from legacy_geometry import dihedral

@njit(cache=True)
def energy(sel,one,pair):
 total=0.
 for i in range(len(sel)):
  total+=one[i,sel[i]]
  for j in range(i):total+=pair[i,sel[i],j,sel[j]]
 return total
@njit(cache=True)
def anneal(nopt,one,pair,seed,steps):
 np.random.seed(seed);sel=np.array([np.random.randint(n) for n in nopt]);v=energy(sel,one,pair);best=v;bs=sel.copy()
 for step in range(steps):
  i=np.random.randint(len(sel));old=sel[i];sel[i]=np.random.randint(nopt[i]);v1=energy(sel,one,pair);temp=3*(.02/3)**(step/max(steps-1,1))
  if v1<v or np.random.random()<math.exp(min(0.,(v-v1)/temp)):
   v=v1
   if v<best:best=v;bs=sel.copy()
  else:sel[i]=old
 sel=bs.copy()
 for sweep in range(4):
  changes=0
  for i in range(len(sel)):
   old=sel[i];best=energy(sel,one,pair);win=old
   for k in range(nopt[i]):
    sel[i]=k;v=energy(sel,one,pair)
    if v<best:best=v;win=k
   sel[i]=win;changes+=old!=win
  if not changes:break
 return sel,energy(sel,one,pair)

def intervals(rows):
 indices=[i for i,r in enumerate(rows) if r['loop']];out=[]
 for i in indices:
  if not out or i!=out[-1][1]:out.append([i,i+1])
  else:out[-1][1]+=1
 assert len(out)==3,out
 return out

def precompute(bb,angles,parent_rows,left,right,target,mouse,gly):
 n=len(bb);allopts=[];fixed=[r for r in parent_rows if not left<r['position']<right]
 fixedcharges=[(CHARGES[r['aa']],charge_center(r['aa'],r['atoms'])) for r in fixed if r['aa'] in 'DERK']
 for i in range(n):
  phi,psi=np.degrees(angles[i]);opts=[];envparts=[]
  for r in fixed:
   names=list(r['atoms'])
   if (i==0 and r['position']==left) or (i==n-1 and r['position']==right):names=[x for x in names if x not in BB]
   envparts.append(atom_data(r['aa'],r['atoms'],names))
  for j in range(n):
   if abs(i-j)>1:envparts.append(atom_data('G',{k:bb[j,z] for z,k in enumerate(BB)}))
  env=np.concatenate(envparts);tree=cKDTree(env[:,:3]);palette='G' if phi>0 else 'AGSTNQDEKVILFYP'
  for aa in palette:
   if aa=='P' and not -95<=phi<=-40:continue
   bestaa=[]
   for rot in range(len(LIB[aa]['xyz'])):
    a=conformer(aa,rot,bb[i]);dat=atom_data(aa,a,LIB[aa]['names'])
    if aa=='P' and not 1.25<=np.linalg.norm(a['N']-a['CD'])<=1.70:continue
    if len(dat) and min(gly.query(dat[:,:3])[0])<2.55:continue
    rp,vw,hb=pair_score(dat,env)
    if rp>18:continue
    targetparts=[pair_score(dat,x) for x in [target,mouse]]
    if max(x[0] for x in targetparts)>10:continue
    q=CHARGES.get(aa,0);center=charge_center(aa,a);elec=0.
    if q and center is not None:
     elec=sum(float(screened_charge(q,q0,np.linalg.norm(center-c0),40.)) for q0,c0 in fixedcharges)
    # Explicitly declared local solvent heuristic, not an SASA/free-energy model.
    solv=0.;contacts=0
    for row in dat:
     close=tree.query_ball_point(row[:3],5.0);k=len(close);contacts+=k
     if row[4]:solv+=.13/(1+.20*k)
     elif k>=8:solv+=.12 # small buried polar penalty; detailed FF follows
    gly_cost=.55 if aa=='G' and phi<0 else 0.
    charge_cost=.16 if aa in 'DEKR' else 0.
    val=rp+vw+1.5*hb+elec+solv+gly_cost+charge_cost+.25*max(sum(x) for x in targetparts)
    bestaa.append({'aa':aa,'rotamer':rot,'atoms':a,'data':dat,'one':float(val),'rep':float(rp),'target_rep':float(max(x[0] for x in targetparts)),'contact_pairs_5A':contacts,'charge_center':center,'q':q})
   # Keep multiple distinct chemistry types; avoid all options being glycine.
   opts.extend(sorted(bestaa,key=lambda x:x['one'])[:2])
  if not opts:return None
  allopts.append(sorted(opts,key=lambda x:x['one']))
 kmax=max(map(len,allopts));one=np.full((n,kmax),1e6);pair=np.zeros((n,kmax,n,kmax));nopt=np.array(list(map(len,allopts)),dtype=np.int64)
 for i,opts in enumerate(allopts):
  for k,o in enumerate(opts):one[i,k]=o['one']
  for j in range(i):
   for a,oi in enumerate(opts):
    for b,oj in enumerate(allopts[j]):
     rep,vw,hb=pair_score(oi['data'],oj['data']);v=rep+vw+1.5*hb
     if oi['q'] and oj['q']:v+=float(screened_charge(oi['q'],oj['q'],np.linalg.norm(oi['charge_center']-oj['charge_center']),40.))
     pair[i,a,j,b]=pair[j,b,i,a]=v
 return allopts,nopt,one,pair

def run(tag,replicates=3):
 t=time.time();out=S5/'intermediate/designs';out.mkdir(parents=True,exist_ok=True);root=S5/'intermediate/backbones';meta=json.loads((root/(tag+'.json')).read_text());bbs=np.load(root/(tag+'.npz'))['bb'];parent,rows,tx=load_parent();left=meta['left_parent_anchor'];right=meta['right_parent_anchor'];n=meta['new_segment_length'];fixedtarget=json.loads((R/'intermediate/human6ARU_aligned.json').read_text());mouse=json.loads((R/'intermediate/mouseAF_aligned.json').read_text());ta=np.concatenate([atom_data(r['aa'],r['atoms']) for r in fixedtarget]);ma=np.concatenate([atom_data(r['aa'],r['atoms']) for r in mouse]);gx=json.loads((R/'stage2/intermediate/target_glycans.json').read_text());gt=cKDTree([a['xyz'] for a in gx]);records=[];seen=set()
 for bidx,bb in enumerate(bbs):
  pc=precompute(bb,np.array(meta['accepted'][bidx]['angles_radians']),rows,left,right,ta,ma,gt)
  if pc is None:records.append({'backbone':bidx,'status':'no compatible rotamer at one or more positions'});continue
  opts,nopt,one,pair=pc
  for rep in range(replicates):
   seed=meta['seed']+bidx*10000+rep*1000+17;sel,value=anneal(nopt,one,pair,seed,8000);chosen=[opts[i][k] for i,k in enumerate(sel)];seqnew=''.join(o['aa'] for o in chosen);new=[]
   for i,o in enumerate(chosen):new.append({'aa':o['aa'],'loop':True,'original_framework_position':None,'parent_position':None,'atoms':{k:np.array(v).tolist() for k,v in o['atoms'].items()},'new_segment_position':i+1})
   full=copy.deepcopy(rows[:left])+new+copy.deepcopy(rows[right-1:]);seq=''.join(r['aa'] for r in full)
   if seq in seen:continue
   seen.add(seq);cid=f'R{2 if meta["region"]=="H2" else 3}{n:02d}{bidx:02d}{rep:02d}';status='candidate'
   if re.search('N[^P][ST]',seq):status='rejected_glycosylation_sequon'
   if len(re.findall(r'[FILVWY]{4}',seq)):status='rejected_hydrophobic_run'
   if max(o['rep'] for o in chosen)>15:status='rejected_packing_repulsion'
   for i,r in enumerate(full,1):r['position']=i
   ints=intervals(full);d=copy.deepcopy(parent);d.update(candidate_id=cid,sequence=seq,structure=full,cdr_intervals_zero_based=ints,cdr_sequences=[seq[a:b] for a,b in ints],campaign_stage=5,branch='shortened_CCD_loop_return',status='unvalidated computational hypothesis',seed=int(seed),historical_rng='MT19937',ancestry={'parent_candidate':'C00003','parent_type':'own untested de novo design','parent_file':str(PARENT),'parent_sha256':sha(PARENT),'relaxed_backbone_source':str(SOURCE),'relaxed_backbone_sha256':sha(SOURCE),'backbone_file':str(root/(tag+'.npz')),'backbone_sha256':sha(root/(tag+'.npz')),'backbone_index':bidx},loop_rebuild={'region':meta['region'],'parent_anchor_positions':[left,right],'parent_deleted_positions':[left+1,right-1],'new_length':n,'old_length':right-left-1,'new_sequence':seqnew,'seed':int(seed),'packing_objective':float(value),'filter_status':status,'method':'CCD2003 + finite empirical rotamers + MT19937 Metropolis annealing; unchanged own-design contact anchors'},metrics={'scores_valid':False,'note':'No inherited affinity/pH ranking score; needs fresh evaluation'})
   for key in ['charge_balance_edits','parent_metrics_not_recomputed','search_corrections']:d.pop(key,None)
   d['sidechain_sources']=[{'position':left+i+1,'aa':o['aa'],'rotamer':o['rotamer'],'source':LIB[o['aa']]['sources'][o['rotamer']]} for i,o in enumerate(chosen)]
   dp=out/(cid+'.json');assert not dp.exists(),dp;dp.write_text(json.dumps(d));records.append({'candidate_id':cid,'backbone_index':bidx,'replicate':rep,'sequence_length':len(seq),'new_sequence':seqnew,'objective':float(value),'status':status,'sequence_sha256':hashlib.sha256(seq.encode()).hexdigest()});print(cid,status,seqnew,round(value,2),flush=True)
 result={'tag':tag,'records':records,'seconds':time.time()-t,'method_qualification':'Search objective only; not folding or binding validation.'};(S5/'reference'/(tag+'_sequence_design.json')).write_text(json.dumps(result,indent=2));print('DONE',tag,len(records),round(result['seconds'],1),flush=True)
if __name__=='__main__':run(sys.argv[1],int(sys.argv[2]) if len(sys.argv)>2 else 3)
