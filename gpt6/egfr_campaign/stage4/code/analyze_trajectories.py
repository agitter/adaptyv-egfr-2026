"""Separate loop internal deformation from motion relative to framework.
Kabsch1976 fits; descriptive time averages, not independent observations or CI.
"""
from pathlib import Path
import json,sys
import numpy as np
R=Path('/mnt/data/egfr_campaign');S=R/'stage4';sys.path.insert(0,str(S/'code'))
from trajectory_contacts import atoms,align
from audit_refined import audit

def root_and_id(cid):
 if cid.startswith('MATCHED_'):return S/'intermediate/matched_controls',cid.removeprefix('MATCHED_')
 return next(r for r in [S,R/'stage3'] if (r/'intermediate/designs'/(cid+'.json')).exists()),cid

def main():
 records=[]
 for p in sorted((S/'intermediate/dynamics').glob('*_v2_s*.json')):
  if p.name.endswith('_progress.json'):continue
  d=json.loads(p.read_text());cid=d['candidate_id'];trace=np.load(p.with_name(p.stem+'_CA_trace.npz'));coords=trace['coords_A'];times=trace['times_ps'];prod=times>=d['warmup_ps'];final=p.with_name(p.stem+'_final.pdb')
  if cid=='NATIVE_3EAK':
   intervals=[[25,38],[52,67],[99,117]];reference=coords[0];qualification='Reference for all internal-deformation comparisons is the initial recorded, minimized unbound frame.';native_sequence=''.join(__import__('Bio.SeqUtils',fromlist=['seq1']).seq1(r.name) for r in __import__('trajectory_contacts').app.PDBFile(str(final)).topology.residues());geom=audit(final,{'sequence':native_sequence,'candidate_id':cid},chain='A')
  else:
   root,actual=root_and_id(cid);design=json.loads((root/'intermediate/designs'/(actual+'.json')).read_text());source=Path(d['metadata']['source']);reference=np.array([r['CA'] for r in atoms(source,'B')]);mask=[r['loop'] for r in design['structure']];intervals=[];a=None
   for i,b in enumerate([*mask,False]):
    if b and a is None:a=i
    elif not b and a is not None:intervals.append([a,i]);a=None
   assert len(intervals)==3;geom=audit(final,design,chain='A');reference=coords[0];qualification='Reference for all internal-deformation comparisons is the initial recorded, minimized unbound frame.'
  vals=[]
  for x in coords:
   row=[]
   for a,b in intervals:
    rr,tt=align(x[a:b],reference[a:b]);row.append(float(np.sqrt(np.mean(np.sum((x[a:b]@rr+tt-reference[a:b])**2,1)))))
   vals.append(row)
  vals=np.array(vals);record={'trajectory':p.stem,'candidate_id':cid,'warmup_ps':d['warmup_ps'],'production_ps':d['production_ps'],'seed':d['seed'],'global_loop_displacement_mean_A':d['summary']['production_mean_loop_RMSD_A'],'loop_internal_deformation_mean_A':vals[prod].mean(0).tolist(),'loop_internal_deformation_final_A':vals[-1].tolist(),'loop_intervals_zero_based_half_open':intervals,'end_point_independent_geometry':geom,'end_point_geometry_qualification':'The omega/disulfide thresholds were designed for minimized structures. Thermal excursions at one instantaneous endpoint are diagnostic only, not automatic exclusion. Chirality inversion or severe bond/clash defects remain separate warnings.','reference_qualification':qualification,'limitations':'Descriptive averages of correlated frames from a short nonequilibrium stress test; not a stability free energy, not statistical significance and not equilibrium folding.'};records.append(record);out=S/'intermediate/trajectory_analysis';out.mkdir(exist_ok=True);(out/(p.stem+'.json')).write_text(json.dumps(record,indent=2));print(p.stem,'internal',np.round(record['loop_internal_deformation_mean_A'],3),'pass',geom['pass'] if geom else None,flush=True)
 (S/'reference/trajectory_analysis_summary.json').write_text(json.dumps(records,indent=2))
if __name__=='__main__':main()
