"""Classical binding-polynomial checks of a conditional proton-linkage bound.
Wyman1964; positive polynomial coefficients. No learned model or pKa predictor.
"""
from pathlib import Path
import json,math
import numpy as np
S=Path('/mnt/data/egfr_campaign/stage4')
def log_partition(coeff,pH):
    x=np.log(np.asarray(coeff,float))-np.arange(len(coeff))*pH*math.log(10)
    m=float(max(x));return m+math.log(float(np.exp(x-m).sum()))
def log_K(bound,free,pH):return log_partition(bound,pH)-log_partition(free,pH)
def main():
    rng=np.random.Generator(np.random.MT19937(19641996));tests=0;max_excess=0.
    for n in range(1,7):
        for k in range(100):
            # Include strongly shifted and nearly neutral protonation equilibria.
            b=np.exp(rng.uniform(-20,120,n+1));f=np.exp(rng.uniform(-20,120,n+1))
            diff=log_K(b,f,6.5)-log_K(b,f,7.4)
            upper=n*.9*math.log(10)
            excess=abs(diff)-upper;max_excess=max(max_excess,excess)
            assert excess<1e-10;tests+=1
            pH=6.9;eps=1e-5
            slope=(log_K(b,f,pH+eps)-log_K(b,f,pH-eps))/(2*eps*math.log(10))
            assert abs(slope)<=n+1e-7;tests+=1
    rt=.00198720425864083*300.
    rows=[{'maximum_linked_protons':n,'maximum_acid_to_neutral_association_constant_ratio':10**(.9*n),'maximum_contrast_kcal_at300K':rt*math.log(10)*.9*n} for n in [1,2,3,4]]
    out={'method':'Positive binding polynomials and proton-linkage identity; all coefficients arbitrary classical equilibrium weights.','test_count':tests,'pass':True,'max_bound_violation':max_excess,'pH_low':6.5,'pH_high':7.4,'temperature_K':300,'bounds':rows,'qualification':'The two-site bound applies ONLY if binding changes proton uptake by at most two protons across this interval. It is not a universal protein limit; other titratable groups, coupled folding and state changes can change the model. The competition detection threshold is unspecified. These are thermodynamic bounds, not predicted candidate affinities.','sources':[{'year':1964,'doi':'10.1016/S0065-3233(08)60190-4'},{'year':1996,'doi':'10.1016/S0006-3495(96)79403-1'}]}
    (S/'reference/proton_linkage_bound_tests.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
if __name__=='__main__':main()
