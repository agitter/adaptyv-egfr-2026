"""Recompute saved CA-trace metrics independently; never use partial trajectories.
Comparisons are descriptive, with two prescribed replicas, not hypothesis tests
or equilibrium folding free energies. All coordinates and source hashes remain.
"""
from pathlib import Path
import sys,json,hashlib
import numpy as np
R=Path('/mnt/data/egfr_campaign');S=R/'stage7'
sys.path[:0]=[str(R/'stage2/code'),str(R/'stage2/runtime/python'),str(R/'code')]
from mm_refine import app,u

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def fit_metrics(x,reference,mask):
 fw=~np.asarray(mask,dtype=bool);c=x[fw].mean(0);c0=reference[fw].mean(0)
 left,_,right=np.linalg.svd((x[fw]-c).T@(reference[fw]-c0));rot=left@right
 if np.linalg.det(rot)<0:left[:,-1]*=-1;rot=left@right
 delta=((x-c)@rot+c0)-reference;d2=np.sum(delta*delta,axis=1)
 boundaries=np.flatnonzero(np.diff(np.r_[False,mask,False]).astype(int))
 return {'loop':float(np.sqrt(d2[~fw].mean())),'framework':float(np.sqrt(d2[fw].mean())),'CDRs':[float(np.sqrt(d2[a:b].mean())) for a,b in zip(boundaries[::2],boundaries[1::2])]}

def one(p):
 d=json.loads(p.read_text());cid=d['candidate_id'];src=Path(d['metadata']['source']);protein=app.PDBFile(str(src));xyz=np.asarray(protein.positions.value_in_unit(u.angstrom))
 chain='A' if cid=='NATIVE_3EAK' else 'B';ids=[a.index for a in protein.topology.atoms() if a.name=='CA' and a.residue.chain.id==chain];reference=xyz[ids]
 if cid=='NATIVE_3EAK':intervals=[[25,38],[52,67],[99,117]];mask=np.array([any(a<=i<b for a,b in intervals) for i in range(len(ids))])
 else:design=json.loads((S/'intermediate/designs'/f'{cid}.json').read_text());mask=np.array([r['loop'] for r in design['structure']])
 trace=p.with_name(p.stem+'_CA_trace.npz');arr=np.load(trace);x=arr['coords_A'];times=arr['times_ps'];frames=d['frames'];assert len(times)==len(frames) and x.shape[1:]==reference.shape and len(mask)==len(reference)
 err=0.;metrics=[]
 for xx,frame,t in zip(x,frames,times):
  assert abs(t-frame['time_ps'])<1e-9;m=fit_metrics(xx,reference,mask);err=max(err,abs(m['loop']-frame['loop_CA_RMSD_A']),abs(m['framework']-frame['framework_CA_RMSD_A']),max(abs(a-b) for a,b in zip(m['CDRs'],frame['CDR_CA_RMSD_A'])));metrics.append(m)
 assert err<1e-8,('Trace metrics do not reproduce stored values',str(p),err)
 prod=np.array([f['phase']=='production' for f in frames]);late=np.array([f['phase']=='production' and f['time_ps']>=d['warmup_ps']+d['production_ps']/2 for f in frames]);early=prod & ~late
 assert abs(times[-1]-(d['warmup_ps']+d['production_ps']))<1e-8 and prod.sum()>1 and d['production_ps']==20
 loop=np.array([m['loop'] for m in metrics]);fw=np.array([m['framework'] for m in metrics]);cdr=np.array([m['CDRs'] for m in metrics]);temp=np.array([f['temperature_K'] for f in frames]);g=d.get('geometry_audit',d.get('final_geometry_audit',{}))
 return {'candidate_id':cid,'seed':d['seed'],'source_result':str(p),'source_sha256':sha(p),'trace_sha256':sha(trace),'coordinate_source_sha256':sha(src),'added_chirality_restraints':d['chirality_restraints'],'position_restraints_production':d['position_restraints_production'],'production_ps':d['production_ps'],'warmup_ps':d['warmup_ps'],'production_frames':int(prod.sum()),'loop_mean_A':float(loop[prod].mean()),'loop_early_mean_A':float(loop[early].mean()),'loop_late_mean_A':float(loop[late].mean()),'loop_final_A':float(loop[-1]),'framework_mean_A':float(fw[prod].mean()),'CDR_mean_A':cdr[prod].mean(0).tolist(),'temperature_mean_K':float(temp[prod].mean()),'temperature_min_K':float(temp[prod].min()),'temperature_max_K':float(temp[prod].max()),'independent_metric_max_difference_A':err,'final_geometry_audit':g,'complete':True}

def main():
 rows=[]
 for p in sorted((S/'intermediate/dynamics').glob('*_v2_s*.json')):
  if p.stem.endswith('_progress'):continue
  rows.append(one(p))
 for p in sorted((S/'intermediate/thermal_unrestrained').glob('*/intermediate/dynamics/*_v2_s*.json')):
  if p.stem.endswith('_progress'):continue
  rows.append(one(p))
 contrasts=[];base={r['seed']:r for r in rows if r['candidate_id']=='M00000' and not r['added_chirality_restraints']}
 for cid in ['M00018','M00021']:
  vals=[r for r in rows if r['candidate_id']==cid and not r['added_chirality_restraints']]
  if len(vals)!=2 or set(base)!={198005,198006}:continue
  pairs=[]
  for r in sorted(vals,key=lambda z:z['seed']):
   b=base[r['seed']];pairs.append({'seed':r['seed'],'loop_mean_reduction_A':b['loop_mean_A']-r['loop_mean_A'],'loop_late_reduction_A':b['loop_late_mean_A']-r['loop_late_mean_A'],'framework_RMSD_change_A':r['framework_mean_A']-b['framework_mean_A'],'temperature_difference_K':r['temperature_mean_K']-b['temperature_mean_K']})
  primary=all(p['loop_mean_reduction_A']>=.15 for p in pairs);late_ok=all(p['loop_late_reduction_A']>=0 for p in pairs);temp_ok=all(abs(p['temperature_difference_K'])<=10 for p in pairs)
  contrasts.append({'candidate_id':cid,'paired_seed_descriptions':pairs,'both_prespecified_mean_reductions_at_least_0p15A':primary,'late_half_not_worse_each_seed':late_ok,'mean_temperatures_within10K':temp_ok,'descriptive_thermal_rule_pass':primary and late_ok and temp_ok,'qualification':'Shared seed labels are protocol matching, not a paired statistical sample of identical atomic noise. Two replicas do not establish significance or equilibrium stability.'})
 out={'replicas':rows,'comparisons_without_extra_restraints':contrasts,'complete_unrestrained_replicas':sum(not r['added_chirality_restraints'] for r in rows),'expected_unrestrained_replicas':8,'all_expected_unrestrained_complete':sum(not r['added_chirality_restraints'] for r in rows)==8,'qualification':'Incomplete trajectories excluded, all complete outcomes retained. No equilibrium folding, melting temperature or expression claim.'}
 (S/'reference/thermal_analysis.json').write_text(json.dumps(out,indent=2));print('Completed replica analyses',len(rows),'unrestrained',out['complete_unrestrained_replicas'])
 for r in rows:print(r['candidate_id'],r['seed'],'added_chirality',r['added_chirality_restraints'],'mean',round(r['loop_mean_A'],4),'late',round(r['loop_late_mean_A'],4),'T',round(r['temperature_mean_K'],2))
 for c in contrasts:print(c['candidate_id'],'thermal_rule',c['descriptive_thermal_rule_pass'],c['paired_seed_descriptions'])
if __name__=='__main__':main()
