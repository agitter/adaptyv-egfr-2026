"""Two-sided histidine recognition: classical directional negative design.

A single acceptor can recognize a neutral histidine tautomer. The additional
term favors simultaneous contacts at both ring nitrogens. It remains a search
heuristic and is tested separately using full three-state proton ensembles.
"""
from gated_design import *
from clamp_geometry import option_features,audit,HIST

@njit(cache=True)
def objective(sel,one,pair,ib,hs,tf,aaidx,positions,clamp,mode):
    base=gated_objective(sel,one,pair,ib,hs,tf,aaidx,positions,0)
    c=np.zeros((clamp.shape[2],2))
    for i in range(len(sel)):
        for h in range(len(c)):
            for n in range(2):c[h,n]=max(c[h,n],clamp[i,sel[i],h,n])
    closure=0.;single=0.;count=0.
    for h in range(len(c)):
        strength=min(c[h,0],c[h,1]);closure+=strength
        single+=max(c[h,0],c[h,1]);count+=min(strength/.4,1.)
    weight=6.+4.*mode
    val=base[0]-weight*closure-.8*single-2.0*count
    return val,base[1],base[2],base[3],base[4],base[5],base[6]

@njit(cache=True)
def search(nopt,one,pair,ib,hs,tf,aaidx,positions,clamp,seed,mode,steps=9000):
    np.random.seed(seed);sel=np.array([np.random.randint(k) for k in nopt]);value=objective(sel,one,pair,ib,hs,tf,aaidx,positions,clamp,mode)[0];best=value;bs=sel.copy()
    for step in range(steps):
        i=np.random.randint(len(sel));old=sel[i];sel[i]=np.random.randint(nopt[i]);nv=objective(sel,one,pair,ib,hs,tf,aaidx,positions,clamp,mode)[0]
        temp=3.2*(.05/3.2)**(step/max(steps-1,1))
        if nv<value or np.random.random()<math.exp(min(0.,(value-nv)/temp)):
            value=nv
            if nv<best:best=nv;bs=sel.copy()
        else:sel[i]=old
    sel=bs.copy()
    for sweep in range(3):
        changes=0
        for i in range(len(sel)):
            old=sel[i];val=objective(sel,one,pair,ib,hs,tf,aaidx,positions,clamp,mode)[0];win=old
            for j in range(nopt[i]):
                sel[i]=j;nv=objective(sel,one,pair,ib,hs,tf,aaidx,positions,clamp,mode)[0]
                if nv<val:val=nv;win=j
            sel[i]=win;changes+=old!=win
        if not changes:break
    return sel,objective(sel,one,pair,ib,hs,tf,aaidx,positions,clamp,mode)

def run(rec,nruns=6):
    start=time.time();pid=rec['pose_id'];out=S3/'intermediate/designs';summ=out/(pid+'_clamp_summary.json')
    if summ.exists():return json.loads(summ.read_text())
    report={'pose_id':pid,'branch':'two_nitrogen_negative_design'}
    try:
        t=candidate_table(rec,False,48,True);clamp=np.zeros((*t['one'].shape,len(HIST),2))
        for i,options in enumerate(t['opts']):
            for j,c in enumerate(options):clamp[i,j]=option_features(c['aa'],c['atoms'])
        maxperN=clamp.max(axis=(0,1));upper=float(np.minimum(maxperN[:,0],maxperN[:,1]).sum())
        report['independent_geometric_clamp_upper_bound']=upper
        unique=set();count=0
        # Even weak geometries remain reported, but expensive restarts are limited.
        if upper<.15:nruns=min(nruns,2)
        for rn in range(nruns):
            seed=20030000+int(pid[1:])*100+rn;mode=rn%3
            sel,metrics=search(t['nopt'],t['one'],t['pair'],t['ib'],t['hs'],t['tf'],t['aaidx'],t['sites'],clamp,seed,mode)
            d=serialize_design(rec,t,sel,metrics,rn,mode);d['candidate_id']=pid+'_C'+str(rn).zfill(2);d['seed']=seed
            if d['sequence'] in unique:continue
            unique.add(d['sequence']);count+=1;d['campaign_stage']=3;d['branch']='two_nitrogen_negative_design'
            d['ancestry']={'pose_record':rec,'CDR_sequence_initialization':'independent MT19937 choices, no donor binding-loop sequence','framework_source':'3EAK generic framework only; original loops removed'}
            d['clamp_geometry']=audit(d['structure']);d['metrics']['ph_energy_contrast_surrogate']=float(metrics[5]);d['metrics']['clamp_search_weight']=6.+4.*mode
            d['search_corrections']=['all gated-search corrections','two-neutral-tautomer bypass counterselection','directional simultaneous ND1/NE2 contacts']
            (out/(d['candidate_id']+'.json')).write_text(json.dumps(d))
        report.update(status='completed',candidate_count=count,site_stats=t['stats'],seconds=time.time()-start)
    except Exception as e:
        import traceback
        report.update(status='failed',error=repr(e),traceback=traceback.format_exc(),seconds=time.time()-start)
    summ.write_text(json.dumps(report,indent=2));print(pid,report['status'],report.get('candidate_count'),round(report.get('independent_geometric_clamp_upper_bound',0),3),round(report['seconds'],1),report.get('error',''),flush=True)
    return report

if __name__=='__main__':
    poses=json.loads((S3/'intermediate/docking/poses.json').read_text());select=json.loads((S3/'reference/clamp_selection.json').read_text());ps={p['pose_id']:p for p in poses}
    for pid in select:run(ps[pid],6)
