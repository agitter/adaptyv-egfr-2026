"""Independent free-pKa assumptions at every enumerated acid and histidine.
Classical binding-polynomial arithmetic, not a pKa predictor or a calibrated
binding free energy. Saves every scenario compactly rather than only extremes.
"""
from pathlib import Path
import sys,json,itertools,math,hashlib
import numpy as np
from scipy.special import logsumexp
R=Path('/mnt/data/egfr_campaign');S=R/'stage5'
sys.path[:0]=[str(R/'stage3/code'),str(R/'stage2/code'),str(R/'stage2/runtime/python'),str(R/'code')]
from proton_ensemble import RT
from acid_histidine_ensemble import polynomial as old_polynomial

def values(states,energies,sites,pkas,fractions,pH):
 states=np.asarray(states);pkas=np.atleast_2d(pkas);fractions=np.atleast_2d(fractions)
 if pkas.shape!=fractions.shape or pkas.shape[1]!=len(sites):raise ValueError('Prior dimension mismatch')
 logs=np.zeros((len(pkas),len(states)))
 for i,site in enumerate(sites):
  state=states[None,:,i];p=pkas[:,i,None];f=fractions[:,i,None]
  if site['kind']=='histidine':
   if np.any((f<=0)|(f>=1)):raise ValueError('Invalid neutral-tautomer fraction')
   logs+=np.where(state==2,math.log(10)*(p-pH),np.where(state==0,np.log(f),np.log1p(-f)))
  elif site['kind']=='carboxylate':logs+=np.where(state>0,math.log(10)*(p-pH)+math.log(.5),0.)
  else:raise ValueError(site['kind'])
 free=logsumexp(logs,axis=1);bound=logsumexp(logs-np.asarray(energies)[None,:]/RT,axis=1)
 return -RT*(bound-free)

def run(path,tag):
 path=Path(path);d=json.loads(path.read_text());sites=d['sites'];states=np.asarray([r['state'] for r in d['microstates']]);choices=[]
 for site in sites:
  if site['kind']=='histidine':choices.append([(p,f) for p in [5.8,6.3,6.8] for f in [.2,.5,.8]])
  else:choices.append([(site['assumed_baseline_free_pKa']+off,.5) for off in [0.,.8,1.6]])
 combinations=np.array(list(itertools.product(*choices)));pkas=combinations[:,:,0];fractions=combinations[:,:,1];results=[];old_error=[]
 for j,physical in enumerate(d['microstates'][0]['energy']['states']):
  energies=np.array([r['energy']['states'][j]['delta_kcal_proxy'] for r in d['microstates']]);lo=values(states,energies,sites,pkas,fractions,6.5);hi=values(states,energies,sites,pkas,fractions,7.4);contrast=hi-lo
  assert np.all(np.abs(contrast)<=len(sites)*.9*RT*math.log(10)+1e-8)
  results.append(np.column_stack([lo,hi,contrast]))
  # Earlier tied-prior scenarios must be an exactly recoverable subset.
  for phis,frac,off in itertools.product([5.8,6.3,6.8],[.2,.5,.8],[0.,.8,1.6]):
   pv=[phis if r['kind']=='histidine' else r['assumed_baseline_free_pKa']+off for r in sites];fv=[frac if r['kind']=='histidine' else .5 for r in sites]
   for ph in [6.5,7.4]:
    old=old_polynomial(d['microstates'],sites,ph,phis,frac,off,j)['conditional_binding_energy_kcal'];new=values(states,energies,sites,pv,fv,ph)[0];old_error.append(abs(old-new))
 assert max(old_error)<1e-9
 results=np.asarray(results);a,b=np.unravel_index(results[:,:,2].argmin(),results[:,:,2].shape);out=S/'intermediate/independent_mixed_priors';out.mkdir(exist_ok=True);file=out/(tag+'.npz');np.savez_compressed(file,site_pKas=pkas,site_HIE_fractions=fractions,energies_acid_neutral_contrast_kcal=results)
 minimum=float(results[:,:,2].min());maximum=float(results[:,:,2].max());record={'source':str(path),'source_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'sites':sites,'microstates':len(states),'scenario_count':int(results.shape[0]*results.shape[1]),'independent_prior_combinations_per_physical_model':len(pkas),'physical_model_parameters':[{k:p[k] for k in ['solute_dielectric','kappa_nm_inverse']} for p in d['microstates'][0]['energy']['states']],'minimum_contrast_kcal':minimum,'maximum_contrast_kcal':maximum,'all_scenarios_acid_on':minimum>0,'worst_case':{'physical_model_index':int(a),'site_pKas':pkas[b].tolist(),'site_HIE_fractions':fractions[b].tolist(),'acid_proxy_kcal':float(results[a,b,0]),'neutral_proxy_kcal':float(results[a,b,1]),'contrast_kcal':float(results[a,b,2])},'previous_tied_prior_maximum_error_kcal':max(old_error),'tied_prior_comparisons':len(old_error),'array_file':str(file),'array_sha256':hashlib.sha256(file.read_bytes()).hexdigest(),'qualification':'Independent assumed free-pKas and tautomers, not measured/predicted pKas. Fixed heavy coordinates and incomplete ionizable-site coverage. Not confidence intervals, independent experiments, calibrated KD or neutral-pH nondetection.'};(out/(tag+'.json')).write_text(json.dumps(record,indent=2));print(tag,'scenarios',record['scenario_count'],'min',minimum,'old_error',max(old_error),flush=True);return record
if __name__=='__main__':run(sys.argv[1],sys.argv[2])
