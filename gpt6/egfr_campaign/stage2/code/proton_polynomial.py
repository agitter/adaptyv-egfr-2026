"""Classical coupled binding polynomials with declared electrostatic approximations.

The observed binding factor is Z_complex/(Z_binder*Z_target), not a difference
between the lowest-energy protonation states. All microstates are enumerated
for interface histidines. Standard proton activity factors appear exactly once.
This is a ranking model, NOT an experimental pH ratio or calibrated KD model.
"""
import math,json,itertools
import numpy as np
from scipy.special import logsumexp
from scipy.spatial import cKDTree
from classical_design import *


def microstate_energies(linear,pairs):
    n=len(linear)
    if n>16:raise ValueError('Too many sites for exhaustive enumeration')
    states=((np.arange(2**n,dtype=np.int64)[:,None]>>np.arange(n))&1).astype(float)
    energies=states@np.array(linear)+.5*np.sum((states@np.array(pairs))*states,axis=1)
    return states,energies


def binding_polynomial(free_linear,free_pairs,bind_linear,bind_pairs,pka,pH,base=0.):
    states,uf=microstate_energies(free_linear,free_pairs)
    ub=uf+states@np.asarray(bind_linear)+.5*np.sum((states@np.asarray(bind_pairs))*states,axis=1)
    activity=math.log(10.)*(states@(np.asarray(pka)-pH))
    lf=activity-uf/RT;lb=activity-ub/RT
    lzf=logsumexp(lf);lzb=logsumexp(lb)
    nf=float(np.exp(lf-lzf)@states.sum(1));nb=float(np.exp(lb-lzb)@states.sum(1))
    return {'binding_free_energy_surrogate':float(base-RT*(lzb-lzf)),'free_protons':nf,'bound_protons':nb,'proton_uptake':nb-nf,'microstates':len(states)}


def extract_charges(rows,partner):
    fixed=[];hist=[]
    for i,r in enumerate(rows):
        c=charge_center(r['aa'],r['atoms'])
        if c is None:continue
        rid=r.get('human_pos',r.get('position',r.get('id',i+1)))
        item={'partner':partner,'position':rid,'xyz':c,'aa':r['aa'],'atoms':r['atoms']}
        if r['aa']=='H':hist.append(item)
        else:item['charge']=CHARGES[r['aa']];fixed.append(item)
    return fixed,hist


def build_model(binder,target,epsilon=20.,free_epsilon=60.,desolv_scale=1.):
    fixed_b,h_b=extract_charges(binder,'B');fixed_t,h_t=extract_charges(target,'T')
    xb=np.array([p for r in binder for p in r['atoms'].values()]);xt=np.array([p for r in target for p in r['atoms'].values()])
    trees={'B':cKDTree(xb),'T':cKDTree(xt)};cloud={'B':xb,'T':xt}
    fixed={'B':fixed_b,'T':fixed_t}
    # Include near-interface histidines, including unbound intrachain coupling.
    hist=[h for h in h_b+h_t if trees['T' if h['partner']=='B' else 'B'].query(h['xyz'])[0]<10.]
    n=len(hist)
    if n>16:raise ValueError(f'{n} interface histidines exceeds exhaustive budget')
    fl=np.zeros(n);bl=np.zeros(n);fp=np.zeros((n,n));bp=np.zeros((n,n));details=[]
    for i,h in enumerate(hist):
        own=h['partner'];other='T' if own=='B' else 'B';pos=h['xyz']
        for a in fixed[own]:
            d=np.linalg.norm(a['xyz']-pos)
            if d<12.:fl[i]+=screened_charge(1.,a['charge'],float(d),free_epsilon)
        cross=0.
        for a in fixed[other]:
            d=np.linalg.norm(a['xyz']-pos)
            if d<18.:cross+=screened_charge(1.,a['charge'],float(d),epsilon)
        # Occluded fraction of the ion hydration surface, using classical
        # point-sampled sphere accessibility. It is an approximate Born penalty.
        count=96;z=1.-2.*(np.arange(count)+.5)/count;phi=np.arange(count)*math.pi*(3.-math.sqrt(5.));r=np.sqrt(1-z*z)
        points=pos+3.0*np.stack([r*np.cos(phi),r*np.sin(phi),z],axis=1)
        d0=trees[own].query(points)[0];d1=trees[other].query(points)[0]
        exposed=d0>2.0;occluded=exposed&(d1<2.8);burial=float(occluded.sum()/max(exposed.sum(),1))
        born=desolv_scale*166.03186*(1./epsilon-1./78.5)*(burial/3.0)*.55
        bl[i]=cross+born
        details.append({'partner':own,'position':h['position'],'intrachain_shift':float(fl[i]),'cross_charge_shift':float(cross),'born_shift':float(born),'hydration_occluded_fraction':burial})
        for j in range(i):
            d=float(np.linalg.norm(pos-hist[j]['xyz']))
            if own==hist[j]['partner']:
                fp[i,j]=fp[j,i]=screened_charge(1.,1.,d,free_epsilon)
            else:bp[i,j]=bp[j,i]=screened_charge(1.,1.,d,epsilon)
    # Base electrostatics: fixed charges on both partners, independent of pH.
    base=0.
    for b in fixed_b:
        for t in fixed_t:
            d=float(np.linalg.norm(b['xyz']-t['xyz']))
            if d<18:base+=screened_charge(b['charge'],t['charge'],d,epsilon)
    return {'free_linear':fl,'free_pairs':fp,'bind_linear':bl,'bind_pairs':bp,'details':details,'base_electrostatic':base,'sites':n,'epsilon':epsilon,'free_epsilon':free_epsilon,'desolv_scale':desolv_scale}


def evaluate_model(model,pka=6.3):
    args=[model[k] for k in ['free_linear','free_pairs','bind_linear','bind_pairs']]
    pk=np.full(model['sites'],pka)
    low=binding_polynomial(*args,pk,6.5,model['base_electrostatic']);high=binding_polynomial(*args,pk,7.4,model['base_electrostatic'])
    contrast=high['binding_free_energy_surrogate']-low['binding_free_energy_surrogate']
    assert abs(contrast)<=model['sites']*RT*math.log(10)*.9+1e-8
    return {'low':low,'high':high,'contrast_kcal_surrogate':contrast,'log10_binding_ratio_surrogate':contrast/(RT*math.log(10.))}


def evaluate_design(binder,target,full=True):
    records=[]
    for eps,desolv in ([(12.,1.),(20.,1.),(32.,1.),(20.,2.)] if full else [(20.,1.)]):
        model=build_model(binder,target,eps,60.,desolv)
        for pka in ([5.8,6.3,6.8] if full else [6.3]):
            res=evaluate_model(model,pka);records.append({'epsilon':eps,'desolv_scale':desolv,'pka':pka,**res})
    base=build_model(binder,target)
    return {'site_count':base['sites'],'site_details':base['details'],'scenarios':records,'worst_contrast_kcal_surrogate':min(r['contrast_kcal_surrogate'] for r in records),'best_contrast_kcal_surrogate':max(r['contrast_kcal_surrogate'] for r in records),'all_scenarios_acid_on':all(r['contrast_kcal_surrogate']>0 for r in records)}


def tests():
    results=[]
    def check(name,condition):
        assert condition,name;results.append(name)
    rng=np.random.Generator(np.random.MT19937(19642003))
    for n in range(1,7):
        z=np.zeros(n);mat=np.zeros((n,n));pk=np.full(n,6.3)
        a=binding_polynomial(z,mat,z,mat,pk,6.5,12.34)
        check('zero linkage '+str(n),abs(a['binding_free_energy_surrogate']-12.34)<1e-10)
        shift=rng.normal(-1.,1.,n)
        a=binding_polynomial(z,mat,shift,mat,pk,6.5)
        ref=float(proton_free_energy(shift,6.5).sum())
        check('independent-site analytic '+str(n),abs(a['binding_free_energy_surrogate']-ref)<1e-10)
        for trial in range(10):
            fl=rng.normal(0,.5,n);bl=rng.normal(-.5,1,n);fp=rng.normal(0,.3,(n,n));fp=(fp+fp.T)/2;np.fill_diagonal(fp,0);bp=fp*.5
            low=binding_polynomial(fl,fp,bl,bp,pk,6.5);high=binding_polynomial(fl,fp,bl,bp,pk,7.4)
            check(f'linked proton bound {n} {trial}',abs(high['binding_free_energy_surrogate']-low['binding_free_energy_surrogate'])<=n*.9*RT*math.log(10)+1e-10)
            delta=.0001
            plus=binding_polynomial(fl,fp,bl,bp,pk,6.5+delta);minus=binding_polynomial(fl,fp,bl,bp,pk,6.5-delta)
            slope=(plus['binding_free_energy_surrogate']-minus['binding_free_energy_surrogate'])/(2*delta)
            check(f'Wyman derivative {n} {trial}',abs(slope-RT*math.log(10)*low['proton_uptake'])<1e-6)
    (S/'reference/proton_polynomial_tests.json').write_text(json.dumps({'passed':len(results),'tests':results},indent=2));print('PASSED',len(results))
if __name__=='__main__':tests()
