"""Classical proton-linkage stress test including designed carboxylates.

Enumerates both neutral-His and neutral-carboxylic-acid tautomers. Fixed heavy
coordinates, explicit Amber99SB charges/OBC. Unbound pKas are declared priors,
not estimated pKas. This is NOT constant-pH dynamics or a calibrated affinity.
"""
from proton_ensemble import *
from fast_obc import parameters

def rotation_between(a,b):
    a=np.asarray(a)/np.linalg.norm(a);b=np.asarray(b)/np.linalg.norm(b);v=np.cross(a,b);c=float(a@b);s=np.linalg.norm(v)
    if s<1e-10:
        if c>0:return np.eye(3)
        raise ValueError('Antiparallel acid oxygen vectors')
    K=np.array([[0.,-v[2],v[1]],[v[2],0.,-v[0]],[-v[1],v[0],0.]])
    return np.eye(3)+K+K@K*((1-c)/(s*s))

def polynomial(microstates,sites,pH,pka_h,tautomer_fraction,acid_offset,model):
    logs=[];energies=[];nh=[]
    for row in microstates:
        log=0.;count=0
        for site,state in zip(sites,row['state']):
            if site['kind']=='histidine':
                if state==2:log+=math.log(10)*(pka_h-pH);count+=1
                else:log+=math.log(tautomer_fraction if state==0 else 1-tautomer_fraction)
            elif state:
                log+=math.log(10)*(site['assumed_baseline_free_pKa']+acid_offset-pH)+math.log(.5);count+=1
        logs.append(log);nh.append(count);energies.append(row['energy']['states'][model]['delta_kcal_proxy'])
    logs=np.asarray(logs);energies=np.asarray(energies);nh=np.asarray(nh);free=logsumexp(logs);bound=logsumexp(logs-energies/RT)
    return {'conditional_binding_energy_kcal':float(-RT*(bound-free)),'proton_uptake':float(np.exp(logs-energies/RT-bound)@nh-np.exp(logs-free)@nh)}

def summarize(records,sites):
    scenarios=[]
    for j,physical in enumerate(records[0]['energy']['states']):
        for phis in [5.8,6.3,6.8]:
            for fraction in [.2,.5,.8]:
                for acid_offset in [0.,.8,1.6]:
                    lo=polynomial(records,sites,6.5,phis,fraction,acid_offset,j);hi=polynomial(records,sites,7.4,phis,fraction,acid_offset,j);contrast=hi['conditional_binding_energy_kcal']-lo['conditional_binding_energy_kcal']
                    assert abs(contrast)<=len(sites)*.9*RT*math.log(10)+1e-8
                    scenarios.append({'solute_dielectric':physical['solute_dielectric'],'kappa_nm_inverse':physical['kappa_nm_inverse'],'assumed_unbound_histidine_pKa':phis,'assumed_unbound_HIE_fraction':fraction,'assumed_acid_pKa_offset':acid_offset,'low_pH':lo,'high_pH':hi,'contrast_kcal':contrast})
    return {'scenarios':scenarios,'minimum_contrast_kcal':min(r['contrast_kcal'] for r in scenarios),'maximum_contrast_kcal':max(r['contrast_kcal'] for r in scenarios),'all_scenarios_acid_on':all(r['contrast_kcal']>0 for r in scenarios),'qualification':'Fixed conformation and declared independent free-pKa priors; not measured pKas, KD or neutral-pH nondetection.'}

def probe_acids(cid,acid_positions=(54,103)):
    source=S3/'intermediate/refined'/(cid+'_human6ARU_relaxed.pdb');out=S3/'intermediate/acid_ensemble';out.mkdir(exist_ok=True);start=time.time();pdb=app.PDBFile(str(source));heavy=app.Modeller(pdb.topology,pdb.positions);heavy.delete([a for a in heavy.topology.atoms() if a.element==app.element.hydrogen]);heavy.topology.createDisulfideBonds(heavy.positions)
    his=identify_sites(heavy);sites=[{**h,'kind':'histidine'} for h in his];acids=[]
    for r in heavy.topology.residues():
        if r.chain.id=='B' and int(r.id) in acid_positions:
            if r.name not in ['ASP','GLU']:raise ValueError('Selected binder site is not a carboxylate')
            acids.append({'residue_index':r.index,'chain':'B','pdb_residue_id':r.id,'kind':'carboxylate','resname':r.name,'assumed_baseline_free_pKa':4.0 if r.name=='ASP' else 4.4})
    assert len(acids)==len(acid_positions);sites+=acids
    assert len(sites)<=6,'Explicit stage5 approval covers at most six sites; never truncate'
    hi={r['residue_index'] for r in his};ac={r['residue_index']:r for r in acids};variants=[None]*heavy.topology.getNumResidues()
    for r in heavy.topology.residues():
        if r.name=='HIS':variants[r.index]='HIP' if r.index in hi else 'HIE'
        if r.index in ac:variants[r.index]='ASH' if r.name=='ASP' else 'GLH'
    random.seed(199020042006);np.random.seed(19902004);heavy.addHydrogens(FF,pH=7.4,variants=variants,platform=mm.Platform.getPlatformByName('CPU'))
    with open(out/(cid+'_all_protonated.pdb'),'w') as f:app.PDBFile.writeFile(heavy.topology,heavy.positions,f)
    base=np.asarray(heavy.positions.value_in_unit(u.nanometer));atom_maps={r.index:{a.name:a.index for a in r.atoms()} for r in heavy.topology.residues()};records=[];charge_checks=[];initial_charge=None
    for state in itertools.product([0,1,2],repeat=len(sites)):
        choice={s['residue_index']:v for s,v in zip(sites,state)};xyz=base.copy();remove=[]
        for site in acids:
            ri=site['residue_index'];v=choice[ri];names=atom_maps[ri];asp=site['resname']=='ASP';on1,on2,cn,hn=('OD1','OD2','CG','HD2') if asp else ('OE1','OE2','CD','HE2')
            if v==2:
                i,j,c,h=map(names.get,[on1,on2,cn,hn]);rot=rotation_between(base[j]-base[c],base[i]-base[c]);xyz[i]=base[j];xyz[j]=base[i];xyz[h]=xyz[j]+rot@(base[h]-base[j])
                assert abs(np.linalg.norm(xyz[h]-xyz[j])-np.linalg.norm(base[h]-base[j]))<1e-10
        mod=app.Modeller(heavy.topology,xyz*u.nanometer)
        for atom in mod.topology.atoms():
            ri=atom.residue.index;v=choice.get(ri)
            if ri in hi and ((v==0 and atom.name=='HD1') or (v==1 and atom.name=='HE2')):remove.append(atom)
            if ri in ac and v==0 and atom.name==('HD2' if ac[ri]['resname']=='ASP' else 'HE2'):remove.append(atom)
        mod.delete(remove);energy=evaluate(mod.topology,mod.positions);records.append({'state':list(state),'energy':energy})
        if state[:len(his)]==tuple([0]*len(his)):
            pars=parameters(mod.topology,mod.positions);q=float(sum(pars[1]));protons=sum(v>0 for v in state[len(his):]);
            if initial_charge is None:initial_charge=q
            delta=q-initial_charge;charge_checks.append({'acid_states':list(state[len(his):]),'charge_increment':delta,'expected':protons,'pass':abs(delta-protons)<1e-5})
        if len(records)%27==0:print(cid,'acid-His states',len(records),'seconds',round(time.time()-start,1),flush=True)
    assert all(r['pass'] for r in charge_checks)
    result={'candidate_id':cid,'source':str(source),'sites':sites,'microstates':records,'microstate_count':len(records),'charge_increment_tests':charge_checks,'proton_model':summarize(records,sites),'seconds':time.time()-start,'method':'Fixed-coordinate classical OBC proton linkage; explicit His and acid tautomers; oxygen type/coordinate exchange for acid tautomers; no equilibrium sampling.'};(out/(cid+'_acid_histidine.json')).write_text(json.dumps(result,indent=2));print(cid,'acid+His contrast',result['proton_model']['minimum_contrast_kcal'],result['proton_model']['maximum_contrast_kcal'],flush=True)

def tests():
    sites=[{'kind':'histidine'},{'kind':'carboxylate','assumed_baseline_free_pKa':4.4}];records=[]
    for s in itertools.product(range(3),repeat=2):records.append({'state':list(s),'energy':{'states':[{'delta_kcal_proxy':-3.,'solute_dielectric':1.,'kappa_nm_inverse':0.}]}})
    tests=[]
    for pH in [6.5,7.4]:
        r=polynomial(records,sites,pH,6.3,.2,.8,0);tests.append({'name':'constant energy '+str(pH),'pass':abs(r['conditional_binding_energy_kcal']+3)<1e-10})
    for row in records:row['energy']['states'][0]['delta_kcal_proxy']=(-4. if row['state'][0]==2 else 0.)+(3. if row['state'][1]>0 else 0.)
    d=1e-4;a=polynomial(records,sites,6.5,6.3,.5,.8,0);b=polynomial(records,sites,6.5+d,6.3,.5,.8,0);c=polynomial(records,sites,6.5-d,6.3,.5,.8,0);tests.append({'name':'Wyman derivative mixed acids and histidines','pass':abs((b['conditional_binding_energy_kcal']-c['conditional_binding_energy_kcal'])/(2*d)-RT*math.log(10)*a['proton_uptake'])<1e-7})
    rot=rotation_between([1.,0.,0.],[-.5,.8660254038,0.]);tests.append({'name':'tautomer rotation orthogonality','pass':float(np.max(abs(rot.T@rot-np.eye(3))))<1e-10});tests.append({'name':'tautomer rotation proper determinant','pass':abs(float(np.linalg.det(rot))-1)<1e-10})
    for t in tests:t['pass']=bool(t['pass'])
    (S3/'reference/acid_histidine_tests.json').write_text(json.dumps(tests,indent=2));assert all(t['pass'] for t in tests);print(len(tests),'acid-His tests passed',flush=True)
if __name__=='__main__':
    if len(sys.argv)>1:probe_acids(sys.argv[1],tuple(map(int,sys.argv[2].split(','))) if len(sys.argv)>2 else (54,103))
    else:tests()
