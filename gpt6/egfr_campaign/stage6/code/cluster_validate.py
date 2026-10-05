"""Test a pre-2011 anchored pair-cluster expansion against archived exact energies.
Historical basis: Zhou et al., PRL 95:148103 (2005); Grigoryan et al.,
PLoS Comput Biol 2:e63 (2006). This is a numerical approximation, not a new
protein generator or a fitted experimental affinity model.
"""
from pathlib import Path
import numpy as np,itertools,json,time,hashlib
R=Path('/mnt/data/egfr_campaign');S=R/'stage6'
def expand(samples,states,order=2):
    states=np.asarray(states,dtype=int);n=states.shape[1];zero=(0,)*n
    coeff={};base=np.asarray(samples[zero]);ans=np.tile(base,(len(states),1))
    for k in range(1,order+1):
        for sites in itertools.combinations(range(n),k):
            for vals in itertools.product([1,2],repeat=k):
                key=tuple(zip(sites,vals));vec=[0]*n
                for i,v in key:vec[i]=v
                value=np.asarray(samples[tuple(vec)])-base
                for j in range(1,k):
                    for subset in itertools.combinations(key,j):value=value-coeff[subset]
                coeff[key]=value
                mask=np.all(states[:,sites]==vals,axis=1);ans[mask]+=value
    return ans

def main():
    paths=sorted((R/'stage5/intermediate').glob('**/*acid_histidine.json'))
    results=[]
    for p in paths:
        d=json.loads(p.read_text());ms=d['microstates'];states=np.asarray([m['state'] for m in ms]);n=states.shape[1]
        if n<5:continue
        exact=np.asarray([[e['delta_kcal_proxy'] for e in m['energy']['states']] for m in ms])
        sample={tuple(st):en for st,en in zip(states,exact) if np.count_nonzero(st)<=2};pred=expand(sample,states);err=pred-exact
        hold=np.count_nonzero(states,axis=1)>2
        res={'source':str(p),'source_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'sites':n,'training_points':len(sample),'held_out_points':int(sum(hold)),'maximum_absolute_energy_error_kcal':float(np.max(abs(err[hold]))),'RMS_energy_error_kcal':float(np.sqrt(np.mean(err[hold]**2))),'maximum_error_by_physical_model_kcal':np.max(abs(err[hold]),axis=0).tolist(),'qualification':'Complete held-out enumeration of this archived model only, not a bound for other structures.'}
        results.append(res);print(p.name,n,len(sample),len(ms),res['maximum_absolute_energy_error_kcal'],flush=True)
    (S/'reference/cluster_validation.json').write_text(json.dumps(results,indent=2))
if __name__=='__main__':main()
