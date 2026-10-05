"""Independent brute-force verification of the inherited finite-state tensor sums."""
from pathlib import Path
import sys,json,itertools,math
import numpy as np
from scipy.special import logsumexp
R=Path('/mnt/data/egfr_campaign');S=R/'stage7';sys.path.insert(0,str(R/'stage6/code'))
from mixed_tensor import energies_at_ph,option_sets,RT

def direct(E,sites,priors,ph,j):
 logs=[];counts=[]
 for states in itertools.product(range(3),repeat=len(sites)):
  lw=0.;nh=0
  for st,site,(pk,f) in zip(states,sites,priors):
   if site['kind']=='histidine':
    if st==2:lw+=math.log(10)*(pk-ph);nh+=1
    else:lw+=math.log(f if st==0 else 1-f)
   elif st!=0:lw+=math.log(.5)+math.log(10)*(pk-ph);nh+=1
  logs.append(lw);counts.append(nh)
 a=np.asarray(logs);cnt=np.asarray(counts);z0=logsumexp(a);zb=logsumexp(a-E[:,j]/RT)
 return -RT*(zb-z0),float(np.exp(a-E[:,j]/RT-zb)@cnt-np.exp(a-z0)@cnt)

def main():
 checks=[]
 for n in [1,2,3,4]:
  sites=[{'kind':'histidine'} if i%2==0 else {'kind':'carboxylate','assumed_baseline_free_pKa':4.4} for i in range(n)]
  ix=np.arange(3**n);E=np.stack([-7+.4*np.sin(ix*.7),-2+np.cos(ix*.3)],axis=1)
  opts=option_sets(sites);choices=[tuple(0 for _ in sites),tuple(len(o)//2 for o in opts),tuple(len(o)-1 for o in opts)]
  for ph in [6.5,7.4]:
   tensor=energies_at_ph(E,sites,ph)
   for ch in choices:
    pri=[o[k] for o,k in zip(opts,ch)]
    for j in range(2):
     g,u=direct(E,sites,pri,ph,j);err=abs(g-tensor[ch+(j,)])
     checks.append({'test':f'brute versus tensor n{n} pH{ph} prior{ch} physical{j}','max_error_kcal':err,'pass':err<1e-10})
     h=1e-4;deriv=(direct(E,sites,pri,ph+h,j)[0]-direct(E,sites,pri,ph-h,j)[0])/(2*h);error=abs(deriv-RT*math.log(10)*u)
     checks.append({'test':f'Wyman derivative n{n} pH{ph} prior{ch} physical{j}','error':error,'pass':error<1e-7})
  const=np.ones_like(E)*-7;lo=energies_at_ph(const,sites,6.5);hi=energies_at_ph(const,sites,7.4)
  checks.append({'test':f'constant energy gives zero pH contrast n{n}','pass':float(np.max(abs(hi-lo)))<1e-10})
  pert=E+5.;base=energies_at_ph(E,sites,6.5);shift=energies_at_ph(pert,sites,6.5)
  checks.append({'test':f'constant offset preserved n{n}','pass':float(np.max(abs(shift-base-5)))<1e-10})
 for c in checks:c['pass']=bool(c['pass'])
 out={'checks':checks,'passed':sum(x['pass'] for x in checks),'count':len(checks),'all_pass':all(x['pass'] for x in checks),'qualification':'Mathematical implementation tests, not validation of the physical energy model.'}
 (S/'reference/polynomial_validation.json').write_text(json.dumps(out,indent=2));print(out['passed'],'/',out['count']);assert out['all_pass']
if __name__=='__main__':main()
