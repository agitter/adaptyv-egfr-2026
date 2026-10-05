"""Two-seed descriptive stress-test analysis with length-matched retained sites.
Frames are correlated observations, never independent replicates. Each run is
20ps, not an equilibrium folding calculation. Endpoint receptor overlays hold
receptor coordinates fixed and cannot establish persistent binding.
"""
from pathlib import Path
import sys,json,hashlib
import numpy as np
R=Path('/mnt/data/egfr_campaign');S=R/'stage5';O=S/'intermediate/dynamics'
sys.path[:0]=[str(R/'stage2/code'),str(R/'stage2/runtime/python'),str(R/'code')]
from mm_refine import app,u
from legacy_geometry import kabsch

IDS=['B00000','N2070600','N2091100','NATIVE_3EAK'];SEEDS=[198005,198006]

def rms(x):return float(np.sqrt(np.mean(np.sum(np.asarray(x)**2,axis=-1))))

def load_reference(cid):
 if cid=='NATIVE_3EAK':
  path=O/'native_heavy.pdb';mask=np.array([any(a<=i<b for a,b in [[25,38],[52,67],[99,117]]) for i in range(127)]);d=None
 else:
  d=json.loads((S/'intermediate/designs'/f'{cid}.json').read_text());path=S/'intermediate/refined'/f'{cid}_human6ARU_relaxed.pdb';mask=np.array([r['loop'] for r in d['structure']])
 p=app.PDBFile(str(path));xyz=np.asarray(p.positions.value_in_unit(u.angstrom));ca=[a.index for a in p.topology.atoms() if a.name=='CA' and a.residue.chain.id=='B'];assert len(ca)==len(mask)
 return xyz[ca],mask,d,p,xyz

def run():
 inputs={cid:load_reference(cid) for cid in IDS if (cid!='NATIVE_3EAK' or (O/'native_heavy.pdb').exists())}
 common=set.intersection(*[{r['parent_position'] for r in inputs[cid][2]['structure'] if r['loop'] and r.get('parent_position') is not None} for cid in IDS if cid!='NATIVE_3EAK'])
 reports=[];missing=[]
 for cid in IDS:
  for seed in SEEDS:
   tag=f'{cid}_v2_s{seed}';path=O/(tag+'.json')
   if not path.exists():missing.append(tag);continue
   meta=json.loads(path.read_text());trace=np.load(O/(tag+'_CA_trace.npz'));coords=trace['coords_A'];times=trace['times_ps'];ref,mask,d,p,xyz=inputs[cid]
   assert coords.shape[1:]==ref.shape and len(meta['frames'])==len(coords)
   framework=~mask;starts=np.flatnonzero(np.diff(np.r_[False,mask,False]));intervals=list(zip(starts[::2],starts[1::2]));production=times>=meta['warmup_ps'];ind=[];site_records=[];loopr=[];fw=[];same=[];agreement=[]
   indexes={r['parent_position']:i for i,r in enumerate(d['structure']) if r.get('parent_position') is not None} if d else {}
   common_indices=[indexes[i] for i in sorted(common)] if d else []
   for t,x in enumerate(coords):
    rot,offset=kabsch(x[framework],ref[framework]);fitted=x@rot+offset;delta=fitted-ref;loopr.append(rms(delta[mask]));fw.append(rms(delta[framework]));agreement.append(abs(loopr[-1]-meta['frames'][t]['loop_CA_RMSD_A']))
    same.append(rms(delta[common_indices]) if d else None)
    local=[]
    for a,b in intervals:
     lr,lo=kabsch(x[a:b],ref[a:b]);local.append(rms(x[a:b]@lr+lo-ref[a:b]))
    ind.append(local);site_records.append({str(i):float(np.linalg.norm(delta[indexes[i]])) for i in [54,100,103,106] if i in indexes})
   assert max(agreement)<2e-9,('Independent alignment disagrees with original record',cid,max(agreement))
   local=np.asarray(ind)[production];result={'candidate_id':cid,'seed':seed,'source':str(path),'source_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'production_ps':meta['production_ps'],'warmup_ps':meta['warmup_ps'],'production_frames':int(production.sum()),'mean_loop_CA_displacement_A':float(np.mean(np.array(loopr)[production])),'mean_framework_CA_displacement_A':float(np.mean(np.array(fw)[production])),'mean_per_loop_internal_deformation_A':np.mean(local,axis=0).tolist(),'mean_retained_loop_CA_displacement_A':float(np.mean(np.array(same,dtype=float)[production])) if d else None,'shared_parent_positions':sorted(common) if d else None,'mean_protected_site_CA_displacement_A':{str(i):float(np.mean([r[str(i)] for j,r in enumerate(site_records) if production[j]])) for i in [54,100,103,106]} if d else {},'mean_temperature_K':meta['summary']['production_mean_temperature_K'],'geometry_pass':meta['final_geometry_audit']['pass'],'maximum_alignment_reimplementation_error_A':max(agreement),'original_summary':meta['summary']}
   # Endpoint side-chain distances after framework alignment to original bound geometry.
   if d:
    endpoint=app.PDBFile(str(O/(tag+'_final.pdb')));ex=np.asarray(endpoint.positions.value_in_unit(u.angstrom));end_ca=[a.index for a in endpoint.topology.atoms() if a.name=='CA'];rot,offset=kabsch(ex[end_ca][framework],ref[framework]);ex=ex@rot+offset;target={int(r.id)+333:{a.name:xyz[a.index] for a in r.atoms()} for r in p.topology.residues() if r.chain.id=='A'};br={int(r.id):{a.name:ex[a.index] for a in r.atoms()} for r in endpoint.topology.residues()};originalbr={int(r.id):{a.name:xyz[a.index] for a in r.atoms()} for r in p.topology.residues() if r.chain.id=='B'};contacts=[]
    for pp,hp in [(54,383),(103,358),(106,358)]:
     i=indexes[pp]+1;oxygen=[a for a in originalbr[i] if a in ['OD1','OD2','OE1','OE2']]
     for hn in ['ND1','NE2']:
      before=min(float(np.linalg.norm(originalbr[i][o]-target[hp][hn])) for o in oxygen);after=min(float(np.linalg.norm(br[i][o]-target[hp][hn])) for o in oxygen);contacts.append({'parent_position':pp,'target_position':hp,'histidine_nitrogen':hn,'initial_min_O_N_A':before,'endpoint_overlay_min_O_N_A':after})
    result['endpoint_fixed_receptor_overlay']=contacts
   reports.append(result)
 aggregates=[]
 for cid in IDS:
  rr=[r for r in reports if r['candidate_id']==cid]
  if not rr:continue
  metrics=['mean_loop_CA_displacement_A','mean_framework_CA_displacement_A','mean_retained_loop_CA_displacement_A','mean_temperature_K'];ag={'candidate_id':cid,'independent_trajectory_count':len(rr),'seeds':[r['seed'] for r in rr]}
  for name in metrics:
   vals=[r[name] for r in rr if r[name] is not None]
   if vals:ag[name]={'mean_of_trajectory_means':float(np.mean(vals)),'range_of_trajectory_means':[min(vals),max(vals)],'values':vals}
  aggregates.append(ag)
 result={'complete':not missing,'missing':missing,'per_trajectory':reports,'aggregates':aggregates,'aggregate_production_ps':sum(r['production_ps'] for r in reports),'aggregate_warmup_plus_production_ps':sum(r['production_ps']+r['warmup_ps'] for r in reports),'qualification':'Two seeds per sequence when complete; descriptive short nonconverged stress tests. No p values from correlated frames; no folding/free-energy/binding inference. Native and designed loop lengths differ. Shared-parent-site comparison does not evaluate the newly rebuilt residues. Endpoint fixed-target overlays are not binding trajectories.'}
 (S/'reference/dynamics_analysis.json').write_text(json.dumps(result,indent=2))
 for a in aggregates:print(a['candidate_id'],a['independent_trajectory_count'],a['mean_loop_CA_displacement_A'],flush=True)
 print('MISSING',missing,flush=True);return result
if __name__=='__main__':run()
