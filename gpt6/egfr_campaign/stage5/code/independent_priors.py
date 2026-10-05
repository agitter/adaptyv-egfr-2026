"""Classical binding-polynomial sensitivity with independent site-specific priors.
Wyman linkage predates 2011. All pKas/tautomer weights here are assumptions,
not predicted or measured quantities. Reuses saved fixed-geometry MM energies.
"""
from pathlib import Path
import sys,json,itertools,math
import numpy as np
from scipy.special import logsumexp
R=Path('/mnt/data/egfr_campaign');S=R/'stage5'
sys.path[:0]=[str(R/'stage3/code'),str(R/'stage2/code'),str(R/'stage2/runtime/python'),str(R/'code')]
from proton_ensemble import conditional_polynomial,RT

def polynomial(states,energies,pkas,fractions,pH):
 states=np.asarray(states);pkas=np.asarray(pkas);fractions=np.asarray(fractions)
 if not (len(pkas)==len(fractions)==states.shape[1]):raise ValueError('Site-prior dimension mismatch')
 if np.any((fractions<=0)|(fractions>=1)):raise ValueError('Tautomer priors must lie in (0,1)')
 logs=np.sum(np.where(states==2,math.log(10)*(pkas-pH),np.where(states==0,np.log(fractions),np.log(1-fractions))),axis=1)
 counts=(states==2).sum(1);lf=logsumexp(logs);lb=logsumexp(logs-energies/RT)
 return {'energy':float(-RT*(lb-lf)),'proton_uptake':float(np.exp(logs-energies/RT-lb)@counts-np.exp(logs-lf)@counts)}

def summarize(source):
 d=json.loads(Path(source).read_text());micro=d['microstates'];states=np.asarray([r['state'] for r in micro]);n=states.shape[1];scenarios=[];uniform_errors=[]
 for j,physical in enumerate(micro[0]['energy']['states']):
  energies=np.array([r['energy']['states'][j]['delta_kcal_proxy'] for r in micro])
  for pkas in itertools.product([5.8,6.3,6.8],repeat=n):
   for frac in itertools.product([.2,.5,.8],repeat=n):
    low=polynomial(states,energies,pkas,frac,6.5);high=polynomial(states,energies,pkas,frac,7.4);contrast=high['energy']-low['energy']
    assert abs(contrast)<=n*.9*RT*math.log(10)+1e-8
    if len(set(pkas))==len(set(frac))==1:
     old=[conditional_polynomial(micro,pkas[0],frac[0],ph,j)['conditional_binding_energy_kcal'] for ph in [6.5,7.4]]
     uniform_errors.extend([abs(old[0]-low['energy']),abs(old[1]-high['energy'])])
    scenarios.append({'solute_dielectric':physical['solute_dielectric'],'kappa_nm_inverse':physical['kappa_nm_inverse'],'assumed_site_pKas':pkas,'assumed_HIE_fractions':frac,'acid_proxy_kcal':low['energy'],'neutral_proxy_kcal':high['energy'],'contrast_kcal':contrast})
 assert max(uniform_errors,default=0.)<1e-10
 result={'source':str(source),'sites':d['sites'],'assumptions':'Independent unbound site priors, not a pKa prediction; conformational entropy omitted','scenario_count':len(scenarios),'minimum_contrast_kcal':min(r['contrast_kcal'] for r in scenarios),'maximum_contrast_kcal':max(r['contrast_kcal'] for r in scenarios),'all_acid_on':all(r['contrast_kcal']>0 for r in scenarios),'uniform_subset_max_error_kcal':max(uniform_errors,default=0.),'scenarios':scenarios}
 out=S/'intermediate/independent_priors';out.mkdir(exist_ok=True);dest=out/Path(source).name;dest.write_text(json.dumps(result,indent=2));return result

def tests():
 rng=np.random.Generator(np.random.MT19937(19641998));checks=[]
 for n in [1,2,3,4]:
  states=np.array(list(itertools.product([0,1,2],repeat=n)));e=rng.normal(size=len(states));p=rng.uniform(5.8,6.8,n);f=rng.uniform(.2,.8,n)
  a=polynomial(states,e,p,f,6.5);b=polynomial(states,e+3.7,p,f,6.5);checks.append(abs(b['energy']-a['energy']-3.7)<1e-10)
  h=1e-4;der=(polynomial(states,e,p,f,6.5+h)['energy']-polynomial(states,e,p,f,6.5-h)['energy'])/(2*h);checks.append(abs(der-RT*math.log(10)*a['proton_uptake'])<1e-7)
  c=polynomial(states,np.full(len(states),-4.2),p,f,6.5);checks.extend([abs(c['energy']+4.2)<1e-10,abs(c['proton_uptake'])<1e-10])
  order=rng.permutation(len(states));d=polynomial(states[order],e[order],p,f,6.5);checks.append(abs(d['energy']-a['energy'])<1e-10)
  swap=np.arange(n)[::-1];d=polynomial(states[:,swap],e,p[swap],f[swap],6.5);checks.append(abs(d['energy']-a['energy'])<1e-10)
  for k in range(30):
   en=rng.normal(0,10,len(states));pa=rng.uniform(4.,9.,n);ta=rng.uniform(.05,.95,n);lo=polynomial(states,en,pa,ta,6.5);hi=polynomial(states,en,pa,ta,7.4);checks.append(abs(hi['energy']-lo['energy'])<=n*.9*RT*math.log(10)+1e-8)
 assert all(checks);r={'checks':len(checks),'passed':sum(checks),'method':'algebraic tests, not biological validation','rng':'MT19937'};(S/'reference/independent_priors_tests.json').write_text(json.dumps(r,indent=2));print(r)

if __name__=='__main__':
 tests()
 for p in sorted((S/'intermediate/refined').glob('*_proton.json')):
  r=summarize(p);print(p.name,r['scenario_count'],round(r['minimum_contrast_kcal'],4),r['all_acid_on'],flush=True)
