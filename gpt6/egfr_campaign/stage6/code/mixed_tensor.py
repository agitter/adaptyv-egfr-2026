"""Exact finite-state proton-binding polynomial by tensor contraction.
This is ordinary finite summation/linear algebra, not an ML model. Free pKas
and tautomer fractions remain explicit assumptions. No conformational entropy
or calibrated Kd is supplied by this calculation.
"""
import numpy as np,itertools,math
RT=.00198720425864083*298.15

def option_sets(sites):
    return [list(itertools.product([5.8,6.3,6.8],[.2,.5,.8])) if s['kind']=='histidine' else [(s['assumed_baseline_free_pKa']+o,.5) for o in [0.,.8,1.6]] for s in sites]

def energies_at_ph(energies,sites,ph):
    E=np.asarray(energies);n=len(sites);assert E.shape[0]==3**n
    shape=tuple(len(x) for x in option_sets(sites));out=np.empty(shape+(E.shape[1],))
    weights=[]
    for site,opts in zip(sites,option_sets(sites)):
        w=[]
        for pk,f in opts:
            r=10.**(pk-ph)
            if site['kind']=='histidine':w.append([f/(1+r),(1-f)/(1+r),r/(1+r)])
            else:w.append([1/(1+r),.5*r/(1+r),.5*r/(1+r)])
        w=np.asarray(w);assert np.max(abs(w.sum(1)-1))<1e-12;weights.append(w)
    for j in range(E.shape[1]):
        baseline=float(E[:,j].min());x=np.exp(-(E[:,j]-baseline)/RT).reshape((3,)*n)
        for w in weights:x=np.tensordot(x,w.T,axes=(0,0))
        assert x.shape==shape and np.all(x>0)
        out[...,j]=baseline-RT*np.log(x)
    return out

def summarize(E,sites):
    lo=energies_at_ph(E,sites,6.5);hi=energies_at_ph(E,sites,7.4);c=hi-lo
    assert np.max(abs(c))<=len(sites)*.9*RT*math.log(10)+1e-8
    ix=np.unravel_index(c.argmin(),c.shape);opts=option_sets(sites)
    return {'minimum_contrast_kcal':float(c.min()),'maximum_contrast_kcal':float(c.max()),'scenario_count':int(c.size),'all_scenarios_acid_on':bool(np.all(c>0)),'worst_case_site_priors':[{'kind':s['kind'],'pKa':opts[i][ix[i]][0],'HIE_fraction_if_histidine':opts[i][ix[i]][1] if s['kind']=='histidine' else None} for i,s in enumerate(sites)],'worst_case_physical_index':int(ix[-1]),'worst_case_acid_proxy_kcal':float(lo[ix]),'worst_case_neutral_proxy_kcal':float(hi[ix]),'qualification':'Exact polynomial for supplied energies and independently varied assumed priors. Not experimental binding or predicted free pKas.'},c
