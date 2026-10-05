"""Enumerate neutral-His tautomers and protonated His in classical MM/GBSA.

Binding polynomials normalize against explicitly declared unbound pKa and
neutral-tautomer priors. They are conditional calculations, NOT measured or
calibrated pKa/KD predictions. All 3^n microstates are retained when n<=5.
The fixed-heavy-coordinate method does not sample conformational entropy.
"""
from pathlib import Path
import sys,json,itertools,math,time,random
import numpy as np
from scipy.special import logsumexp
R=Path('/mnt/data/egfr_campaign');S3=R/'stage3'
sys.path[:0]=[str(S3/'code'),str(R/'stage2/code'),str(R/'stage2/runtime/python'),str(R/'code')]
from mm_refine import mm,app,u,FF
from fast_obc import evaluate
RT=.0019872041*298.15


def identify_sites(heavy,target_chains=('A',),distance_A=8.):
    xyz=np.asarray(heavy.positions.value_in_unit(u.angstrom));atoms=list(heavy.topology.atoms());part=np.array([a.residue.chain.id not in target_chains for a in atoms]);sites=[]
    for res in heavy.topology.residues():
        if res.name not in ['HIS','HIE','HID','HIP']:continue
        names={a.name:a.index for a in res.atoms()}
        if not {'ND1','NE2'}<=set(names):continue
        center=xyz[[names['ND1'],names['NE2']]].mean(0);own=res.chain.id not in target_chains
        distance=float(np.linalg.norm(xyz[part!=own]-center,axis=1).min())
        if distance<distance_A:sites.append({'residue_index':res.index,'chain':res.chain.id,'pdb_residue_id':res.id,'partner':'binder' if own else 'target','nearest_partner_A':distance})
    if len(sites)>5:raise ValueError('More than five interface histidines; explicit larger-state handling required, not silently truncated')
    return sites


def conditional_polynomial(microstates,pka,tautomer_hie_fraction,pH,model_index):
    n=len(microstates[0]['state']);logs=[];energies=[];counts=[]
    for row in microstates:
        log=0.;nh=0
        for s in row['state']:
            if s==2:log+=math.log(10)*(pka-pH);nh+=1
            else:log+=math.log(tautomer_hie_fraction if s==0 else 1-tautomer_hie_fraction)
        logs.append(log);energies.append(row['energy']['states'][model_index]['delta_kcal_proxy']);counts.append(nh)
    logs=np.array(logs);energies=np.array(energies);counts=np.array(counts)
    zf=logsumexp(logs);zb=logsumexp(logs-energies/RT)
    free=float(np.exp(logs-zf)@counts);bound=float(np.exp(logs-energies/RT-zb)@counts)
    return {'conditional_binding_energy_kcal':float(-RT*(zb-zf)),'free_H_count':free,'bound_H_count':bound,'proton_uptake':bound-free}


def summarize_microstates(microstates):
    scenarios=[];n=len(microstates[0]['state'])
    for j,physical in enumerate(microstates[0]['energy']['states']):
        for pka in [5.8,6.3,6.8]:
            for frac in [.2,.5,.8]:
                lo=conditional_polynomial(microstates,pka,frac,6.5,j);hi=conditional_polynomial(microstates,pka,frac,7.4,j)
                contrast=hi['conditional_binding_energy_kcal']-lo['conditional_binding_energy_kcal']
                if abs(contrast)>n*.9*RT*math.log(10)+1e-8:raise AssertionError('Proton-linkage thermodynamic bound violated')
                scenarios.append({'solute_dielectric':physical['solute_dielectric'],'kappa_nm_inverse':physical['kappa_nm_inverse'],'assumed_unbound_pKa':pka,'assumed_unbound_HIE_fraction':frac,'low_pH':lo,'high_pH':hi,'contrast_kcal':contrast,'conditional_log10_affinity_ratio':contrast/(RT*math.log(10))})
    return {'scenarios':scenarios,'minimum_contrast_kcal':min(s['contrast_kcal'] for s in scenarios),'maximum_contrast_kcal':max(s['contrast_kcal'] for s in scenarios),'all_scenarios_acid_on':all(s['contrast_kcal']>0 for s in scenarios),'maximum_contrast_from_proton_count_kcal':n*.9*RT*math.log(10),'caveat':'Unbound pKas/tautomer populations are assumed priors; neither absolute affinity nor neutral-pH nondetection is established.'}


def probe(pdb_path,out_prefix,target_chains=('A',),target_map=None,forced_sites=None):
    start=time.time();pdb=app.PDBFile(str(pdb_path));heavy=app.Modeller(pdb.topology,pdb.positions)
    heavy.delete([a for a in heavy.topology.atoms() if a.element==app.element.hydrogen]);heavy.topology.createDisulfideBonds(heavy.positions)
    sites=identify_sites(heavy,target_chains)
    if forced_sites is not None:
        required={(str(c),str(i)) for c,i in forced_sites};chosen=[]
        xyz=np.asarray(heavy.positions.value_in_unit(u.angstrom));allatoms=list(heavy.topology.atoms())
        for res in heavy.topology.residues():
            if (res.chain.id,res.id) not in required:continue
            if res.name not in ['HIS','HIE','HID','HIP']:raise ValueError('Forced site is not histidine')
            names={a.name:a.index for a in res.atoms()};center=xyz[[names['ND1'],names['NE2']]].mean(0);own=res.chain.id not in target_chains;other=[a.index for a in allatoms if (a.residue.chain.id not in target_chains)!=own]
            chosen.append({'residue_index':res.index,'chain':res.chain.id,'pdb_residue_id':res.id,'partner':'binder' if own else 'target','nearest_partner_A':float(np.linalg.norm(xyz[other]-center,axis=1).min())})
        if {(r['chain'],r['pdb_residue_id']) for r in chosen}!=required:raise ValueError('Some forced sites not found')
        if len(chosen)>5:raise ValueError('More than five forced sites not supported')
        sites=chosen
    indices={h['residue_index'] for h in sites}
    for h in sites:
        if target_map and h['partner']=='target':h['human_position']=target_map.get(int(h['pdb_residue_id']))
    variants=[None]*heavy.topology.getNumResidues()
    for res in heavy.topology.residues():
        if res.name in ['HIS','HIE','HID','HIP']:variants[res.index]='HIP' if res.index in indices else 'HIE'
    random.seed(19641998);np.random.seed(19641998);heavy.addHydrogens(FF,pH=7.4,variants=variants,platform=mm.Platform.getPlatformByName('CPU'))
    out_prefix=Path(out_prefix);out_prefix.parent.mkdir(parents=True,exist_ok=True)
    with open(str(out_prefix)+'_allHIP.pdb','w') as f:app.PDBFile.writeFile(heavy.topology,heavy.positions,f)
    microstates=[];n=len(sites)
    for state in itertools.product([0,1,2],repeat=n):
        choice={s['residue_index']:v for s,v in zip(sites,state)};mod=app.Modeller(heavy.topology,heavy.positions)
        remove=[]
        for atom in mod.topology.atoms():
            v=choice.get(atom.residue.index,2)
            if (v==0 and atom.name=='HD1') or (v==1 and atom.name=='HE2'):remove.append(atom)
        mod.delete(remove)
        energy=evaluate(mod.topology,mod.positions,target_chains)
        microstates.append({'state':list(state),'variants':['HIE','HID','HIP'],'energy':energy})
    result={'source':str(pdb_path),'site_selection':'explicit common site set' if forced_sites is not None else '8-Angstrom interface cutoff','forced_sites':forced_sites,'sites':sites,'microstate_count':len(microstates),'microstates':microstates,'proton_model':summarize_microstates(microstates),'seconds':time.time()-start,'method':'Amber99SB/OBC2, matched fixed-coordinate subtraction; explicit 3-state histidines; dielectric/salt/unbound-pKa/tautomer sensitivity'}
    Path(str(out_prefix)+'_proton.json').write_text(json.dumps(result,indent=2));print(out_prefix.name,'states',len(microstates),'pH contrast',result['proton_model']['minimum_contrast_kcal'],result['proton_model']['maximum_contrast_kcal'],'seconds',round(time.time()-start,1),flush=True)
    return result


def tests():
    records=[]
    for n in range(1,5):
        states=[]
        for state in itertools.product([0,1,2],repeat=n):states.append({'state':list(state),'energy':{'states':[{'delta_kcal_proxy':-7.,'solute_dielectric':1.,'kappa_nm_inverse':0.}]}})
        for pH in [6.5,7.4]:
            r=conditional_polynomial(states,6.3,.3,pH,0);records.append({'test':f'constant energy n={n} pH={pH}','pass':abs(r['conditional_binding_energy_kcal']+7.)<1e-10})
        s=summarize_microstates(states);records.append({'test':f'zero linkage n={n}','pass':abs(s['minimum_contrast_kcal'])<1e-10})
        for row in states:row['energy']['states'][0]['delta_kcal_proxy']=-2.*row['state'].count(2)
        s=summarize_microstates(states);records.append({'test':f'acid on n={n}','pass':s['minimum_contrast_kcal']>0})
        delta=1e-4;a=conditional_polynomial(states,6.3,.3,6.5,0);b=conditional_polynomial(states,6.3,.3,6.5+delta,0);c=conditional_polynomial(states,6.3,.3,6.5-delta,0)
        derivative=(b['conditional_binding_energy_kcal']-c['conditional_binding_energy_kcal'])/(2*delta)
        records.append({'test':f'Wyman derivative n={n}','pass':abs(derivative-RT*math.log(10)*a['proton_uptake'])<1e-7})
    out={'tests':records,'passed':sum(r['pass'] for r in records),'all_pass':all(r['pass'] for r in records)};(S3/'reference/proton_ensemble_tests.json').write_text(json.dumps(out,indent=2));assert out['all_pass'];print(out['passed'],'proton-ensemble tests passed')

if __name__=='__main__':
    if len(sys.argv)==1:tests()
    else:probe(sys.argv[1],sys.argv[2])
