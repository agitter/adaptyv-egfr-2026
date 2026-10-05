"""Third-order classical cluster expansion of explicit acid/His microstate energies.
All constant, point, pair and triplet basis states are evaluated with Amber99SB/
OBC, then compared with independent held-out microstates. Exact finite-sum pH
polynomials use the interpolated energy tensor. Not exhaustive energy evaluation,
not a learned protein generator, not a calibrated binding free energy.

Uses the same topology/tautomer construction as the archived exhaustive stage5
acid_histidine_six_sites.py. CE precedent: PRL95:148103(2005), PLoSCB2:e63(2006).
"""
from pathlib import Path
import sys,json,itertools,time,hashlib,random
import numpy as np
R=Path('/mnt/data/egfr_campaign');S=R/'stage6'
sys.path[:0]=[str(S/'code'),str(R/'stage5/code'),str(R/'stage3/code'),str(R/'stage2/code'),str(R/'code'),str(R/'stage2/runtime/python')]
from acid_histidine_six_sites import rotation_between
from proton_ensemble import app,mm,u,FF,identify_sites,evaluate
from fast_obc import parameters
from cluster_validate import expand
from mixed_tensor import summarize

def main(cid):
 start=time.time();source=S/'intermediate/refined'/f'{cid}_human6ARU_relaxed.pdb';out=S/'intermediate/cluster_acids';out.mkdir(exist_ok=True)
 final=out/f'{cid}_cluster.json'
 if final.exists():print('CACHED',cid);return
 pdb=app.PDBFile(str(source));heavy=app.Modeller(pdb.topology,pdb.positions);heavy.delete([a for a in heavy.topology.atoms() if a.element==app.element.hydrogen]);heavy.topology.createDisulfideBonds(heavy.positions)
 his=identify_sites(heavy);sites=[{**h,'kind':'histidine'} for h in his];acids=[]
 acid_positions={52,53,54,103,104,105,106}
 for r in heavy.topology.residues():
  if r.chain.id=='B' and int(r.id) in acid_positions and r.name in ['ASP','GLU']:
   acids.append({'residue_index':r.index,'chain':'B','pdb_residue_id':r.id,'kind':'carboxylate','resname':r.name,'assumed_baseline_free_pKa':4.0 if r.name=='ASP' else 4.4})
 sites+=acids;n=len(sites);assert 4<=n<=8,(cid,n)
 hi={r['residue_index'] for r in his};ac={r['residue_index']:r for r in acids};variants=[None]*heavy.topology.getNumResidues()
 for r in heavy.topology.residues():
  if r.name=='HIS':variants[r.index]='HIP' if r.index in hi else 'HIE'
  if r.index in ac:variants[r.index]='ASH' if r.name=='ASP' else 'GLH'
 random.seed(199020042006);np.random.seed(19902004);heavy.addHydrogens(FF,pH=7.4,variants=variants,platform=mm.Platform.getPlatformByName('CPU'))
 with open(out/f'{cid}_all_protonated.pdb','w') as f:app.PDBFile.writeFile(heavy.topology,heavy.positions,f)
 base=np.asarray(heavy.positions.value_in_unit(u.nanometer));atom_maps={r.index:{a.name:a.index for a in r.atoms()} for r in heavy.topology.residues()};evaluated={};qzero=None
 journal=out/f'{cid}_exact_states.jsonl'
 if journal.exists():
  for line in journal.read_text().splitlines():
   z=json.loads(line);evaluated[tuple(z['state'])]=z
 def energy(state,role):
  nonlocal qzero
  if tuple(state) in evaluated:return evaluated[tuple(state)]
  choice={s['residue_index']:v for s,v in zip(sites,state)};xyz=base.copy();remove=[]
  for site in acids:
   ri=site['residue_index'];v=choice[ri];names=atom_maps[ri];asp=site['resname']=='ASP';on1,on2,cn,hn=('OD1','OD2','CG','HD2') if asp else ('OE1','OE2','CD','HE2')
   if v==2:
    i,j,c,h=map(names.get,[on1,on2,cn,hn]);rot=rotation_between(base[j]-base[c],base[i]-base[c]);xyz[i]=base[j];xyz[j]=base[i];xyz[h]=xyz[j]+rot@(base[h]-base[j]);assert abs(np.linalg.norm(xyz[h]-xyz[j])-np.linalg.norm(base[h]-base[j]))<1e-10
  mod=app.Modeller(heavy.topology,xyz*u.nanometer)
  for atom in mod.topology.atoms():
   ri=atom.residue.index;v=choice.get(ri)
   if ri in hi and ((v==0 and atom.name=='HD1') or (v==1 and atom.name=='HE2')):remove.append(atom)
   if ri in ac and v==0 and atom.name==('HD2' if ac[ri]['resname']=='ASP' else 'HE2'):remove.append(atom)
  mod.delete(remove);en=evaluate(mod.topology,mod.positions);row={'state':list(state),'energy':en,'evaluation_role':role}
  if sum(v!=0 for v in state)<=1:
   q=float(sum(parameters(mod.topology,mod.positions)[1]));protons=sum(v==2 for v in state[:len(his)])+sum(v>0 for v in state[len(his):]);row.update(net_charge=q,added_protons=protons)
  evaluated[tuple(state)]=row
  with journal.open('a') as f:f.write(json.dumps(row)+'\n')
  if len(evaluated)%50==0:print(cid,'exact states',len(evaluated),'s',round(time.time()-start,1),flush=True)
  return row
 states=np.asarray(list(itertools.product(range(3),repeat=n)),int);basis=states[np.count_nonzero(states,axis=1)<=3]
 for state in basis:energy(tuple(map(int,state)),'cluster_basis_order3')
 sample={tuple(st):np.asarray([e['delta_kcal_proxy'] for e in evaluated[tuple(st)]['energy']['states']]) for st in basis};pred=expand(sample,states,3)
 # Held-out points are not used in fitting. Select before evaluating their energies.
 rest=np.flatnonzero(np.count_nonzero(states,axis=1)>3);rng=np.random.RandomState(20060000+int(cid[1:]));hold=set(map(int,rng.choice(rest,min(64,len(rest)),replace=False)))
 # Include every all-protonated tautomer combination, beyond the fit support.
 for idx,st in enumerate(states):
  if np.all(st[:len(his)]==2) and np.all(st[len(his):]>0):hold.add(idx)
 # Include extreme predicted-energy states as an additional non-random stress test.
 for j in range(pred.shape[1]):
  for idx in np.argsort(pred[:,j])[:5]:
   if np.count_nonzero(states[idx])>3:hold.add(int(idx))
 errors=[]
 for ix in sorted(hold):
  st=tuple(map(int,states[ix]));exact=np.array([x['delta_kcal_proxy'] for x in energy(st,'held_out_validation')['energy']['states']]);errors.append({'state':list(st),'prediction_minus_exact_kcal':(pred[ix]-exact).tolist()})
 charges=[v for v in evaluated.values() if 'net_charge' in v];zero=next(v['net_charge'] for v in charges if not any(v['state']));assert all(abs(v['net_charge']-zero-v['added_protons'])<1e-5 for v in charges)
 summary,contrasts=summarize(pred,sites);er=np.asarray([r['prediction_minus_exact_kcal'] for r in errors]);mx=float(abs(er).max());np.savez_compressed(out/f'{cid}_energy_tensor.npz',states=states,energies_kcal=pred,contrast_kcal=contrasts)
 result={'candidate_id':cid,'source':str(source),'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'sites':sites,'state_count':len(states),'cluster_order':3,'basis_states':len(basis),'held_out_states':len(hold),'exact_energy_evaluations':len(evaluated),'held_out_max_absolute_energy_error_kcal':mx,'held_out_RMS_energy_error_kcal':float(np.sqrt(np.mean(er**2))),'validation_errors':errors,'charge_increment_tests':len(charges),'charge_increment_pass':True,'validation_threshold_kcal':.05,'numerical_validation_pass':mx<=.05,'proton_model':summary,'seconds':time.time()-start,'qualifications':['Third-order cluster interpolation, not exhaustive microstate energy evaluation.','Held-out errors are observed errors, not a rigorous bound on every untested state.','All acidic residues at positions52,53,54,103,104,105,106 are included when present; other ionizable groups remain outside this scope.','Independent assumed unbound pKas and tautomers; fixed heavy geometry; no affinity or neutral-pH nondetection claim.']}
 final.write_text(json.dumps(result,indent=2));print('FINISHED',cid,'contrast',summary['minimum_contrast_kcal'],'validation',mx,'seconds',round(time.time()-start,1),flush=True)
if __name__=='__main__':main(sys.argv[1])
