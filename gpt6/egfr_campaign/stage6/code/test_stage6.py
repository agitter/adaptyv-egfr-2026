"""Numerical tests. Passing them does not validate biological predictions."""
from pathlib import Path
import numpy as np,json,itertools,math,sys
from scipy.special import logsumexp
from mixed_tensor import energies_at_ph,summarize,option_sets,RT
from cluster_validate import expand
S=Path('/mnt/data/egfr_campaign/stage6');R=S.parent;results=[]
def check(name,ok,**data):
 results.append({'name':name,'pass':bool(ok),**data})
 if not ok:raise AssertionError((name,data))
rng=np.random.RandomState(20061964)
# Tensor contraction versus an independent direct log-sum-exp summation.
for n in range(1,6):
 sites=[{'kind':'histidine'} if i%2==0 else {'kind':'carboxylate','assumed_baseline_free_pKa':4.4} for i in range(n)];states=np.asarray(list(itertools.product(range(3),repeat=n)));E=rng.normal(size=(len(states),2))*4.-30
 for ph in [6.5,7.4]:
  ans=energies_at_ph(E,sites,ph);opts=option_sets(sites)
  for trial in range(12):
   ix=tuple(int(rng.randint(len(o))) for o in opts);lw=np.zeros(len(states))
   for i,(site,op) in enumerate(zip(sites,opts)):
    pk,f=op[ix[i]]
    if site['kind']=='histidine':lw+=np.where(states[:,i]==2,math.log(10)*(pk-ph),np.where(states[:,i]==0,math.log(f),math.log(1-f)))
    else:lw+=np.where(states[:,i]>0,math.log(10)*(pk-ph)+math.log(.5),0.)
   exact=-RT*(logsumexp(lw[:,None]-E/RT,axis=0)-logsumexp(lw));error=float(np.max(abs(exact-ans[ix])));check(f'tensor_direct_n{n}_ph{ph}_{trial}',error<1e-10,error=error)
  flat=np.full_like(E,-7.5);v=energies_at_ph(flat,sites,ph);check(f'flat_energy_n{n}_ph{ph}',np.max(abs(v+7.5))<1e-10)
 summary,c=summarize(E,sites);check(f'linkage_bound_n{n}',np.max(abs(c))<=n*.9*RT*math.log(10)+1e-8)
# Exact recovery of a known order-three function, including signs/cross terms.
for n in [4,5,6]:
 states=np.asarray(list(itertools.product(range(3),repeat=n)));E=np.tile(np.array([2.,-3.]),(len(states),1))
 for k in [1,2,3]:
  for idx in itertools.combinations(range(n),k):
   for vals in itertools.product([1,2],repeat=k):E[np.all(states[:,idx]==vals,axis=1)]+=rng.normal(size=2)
 sample={tuple(st):e for st,e in zip(states,E) if np.count_nonzero(st)<=3};pred=expand(sample,states,3);err=float(np.max(abs(E-pred)));check(f'cluster_order3_n{n}',err<1e-10,error=err)
 # A true four-body term must not be falsely called exact by the third-order fit.
 E[np.all(states[:,:4]==2,axis=1)]+=1.;pred=expand({tuple(st):e for st,e in zip(states,E) if np.count_nonzero(st)<=3},states,3);check(f'cluster_detects_missing_order4_n{n}',np.max(abs(E-pred))>.99)
# Recover archived mixed-prior minima, allowing only the tiny gas-constant
# decimal precision difference relative to the older implementation.
for d in json.loads((S/'reference/cluster_validation_expanded.json').read_text()):
 p=Path(d['source']);old=json.loads(p.read_text());E=np.asarray([[e['delta_kcal_proxy'] for e in m['energy']['states']] for m in old['microstates']]);ans,C=summarize(E,old['sites']);check('archived_exact_recovery_'+str(p),abs(ans['minimum_contrast_kcal']-d['exact_summary']['minimum_contrast_kcal'])<1e-10)
(S/'reference/stage6_tests.json').write_text(json.dumps({'count':len(results),'passed':sum(r['pass'] for r in results),'tests':results,'qualification':'Implementation tests only; no biological success inferred.'},indent=2));print('TESTS',len(results),'PASSED')
