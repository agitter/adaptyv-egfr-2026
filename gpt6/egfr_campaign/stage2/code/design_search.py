"""Build discrete energy tables and optimize CDRs from scratch for each pose."""
from pathlib import Path
import sys,os,json,time,math,hashlib
import numpy as np
from scipy.spatial import cKDTree
from classical_design import *

PRIOR={'A':.9,'C':10.,'D':.65,'E':.8,'F':1.6,'G':1.25,'H':1.,'I':1.4,'K':1.3,'L':1.3,'M':2.,'N':.65,'P':1.0,'Q':.75,'R':1.4,'S':.35,'T':.6,'V':1.2,'W':1.6,'Y':.75}

def target_context():
    targets=[json.loads((O/(name+'_aligned.json')).read_text()) for name in ['human6ARU','mouseAF']]
    hs=sorted(set(r['human_pos'] for rr in targets for r in rr if r['aa']=='H' and 330<=r['human_pos']<=515))
    info=[]
    for rows in targets:
        atoms=all_atoms(rows);groups=[];centers=[];charges=[];hcent={}
        for r in rows:
            c=charge_center(r['aa'],r['atoms'])
            if c is None:continue
            if r['aa']=='H' and r['human_pos'] in hs:hcent[r['human_pos']]=c
            elif r['aa']!='H':centers.append(c);charges.append(CHARGES[r['aa']])
        info.append({'atoms':atoms,'tree':cKDTree(atoms[:,:3]),'centers':np.array(centers),'charges':np.array(charges),'hist':hcent})
    return targets,hs,info

def candidate_table(rec,maxoptions=42):
    bb,mask,seq,orig,pp,fixed=prepare_pose(rec)
    targets,hpos,targ=target_context();sites=np.flatnonzero(mask)
    fixrows=[{'aa':seq[i],'atoms':a} for i,a in fixed.items()]
    allfixed=all_atoms(fixrows)
    opts=[];stats=[]
    for site in sites:
        env=[]
        for j,a in fixed.items():
            if np.linalg.norm(bb[site,1]-bb[j,1])<15:env.append(atom_data(seq[j],a))
        for j in sites:
            if j!=site and np.linalg.norm(bb[site,1]-bb[j,1])<15:env.append(atom_data('A',{n:p for n,p in zip(BB,bb[j])}))
        env=np.concatenate(env) if env else np.zeros((0,10))
        local=[t['atoms'][t['tree'].query_ball_point(bb[site,1],14.)] for t in targ]
        phi,psi=np.degrees(pp[site]);allowed='ADEFGHIKLNPQRSTVWY'
        if phi>20:allowed='G' # non-glycine positive-phi escape is disallowed here
        elif phi<-145 or psi>175:allowed='ADGNSYHETQ'
        cache=[]
        for aa in allowed:
            if aa=='P' and not (-90<phi<-35):continue
            for rot in range(len(LIB[aa]['xyz'])):
                atoms=conformer(aa,rot,bb[site]);dat=atom_data(aa,atoms,LIB[aa]['names'])
                rep,vdw,hb=pair_score(dat,env)
                if rep>18.:continue
                one=rep+vdw+hb+PRIOR[aa]
                if aa=='G' and phi>20:one-=1.0
                if len(dat):
                    ds=np.linalg.norm(dat[:,None,:3]-env[None,:,:3],axis=2)
                    exposed=(ds<5.).sum(1)<4
                    one+=.10*float(np.sum(exposed*(dat[:,4]>0)))
                    # Unsatisfied buried heteroatoms have a modest regularizer.
                    buried=(ds<4.2).sum(1)>6
                    one+=.35*float(np.sum(buried*(dat[:,4]==0)))
                ib=[];hshift=[];fields=[];sterics=[]
                center=charge_center(aa,atoms);q=CHARGES.get(aa,0.)
                for t,envt in zip(targ,local):
                    rp,vw,hb=pair_score(dat,envt);sterics.append(rp)
                    value=rp+vw+hb
                    hs=0.;tf=np.zeros(len(hpos))
                    if center is not None:
                        dist=np.linalg.norm(t['centers']-center,axis=1)
                        potential=sum(screened_charge(1.,float(c),float(d)) for c,d in zip(t['charges'],dist) if d<18)
                        if aa=='H':
                            burial=float(np.sum(np.linalg.norm(envt[:,:3]-center,axis=1)<5.0))
                            hs=potential+min(1.8,.12*burial)
                        else:value+=q*potential
                        if aa!='H':
                            for k,h in enumerate(hpos):
                                if h in t['hist']:
                                    d=np.linalg.norm(center-t['hist'][h]);tf[k]+=screened_charge(q,1.,float(d)) if d<18 else 0.
                    if len(dat):
                        for k,h in enumerate(hpos):
                            if h in t['hist']:
                                d=np.linalg.norm(dat[:,:3]-t['hist'][h],axis=1);tf[k]+=.05*float(np.sum(d<5.0))
                    ib.append(value);hshift.append(hs);fields.append(tf)
                if max(sterics)>22.:continue
                acid=np.array(ib)+proton_free_energy(hshift,6.5)*(aa=='H')+proton_free_energy(np.array(fields),6.5).sum(1)
                high=ib[0]+float(proton_free_energy(hshift[0],7.4))*(aa=='H')+float(proton_free_energy(np.array(fields)[0],7.4).sum())
                pre=one+.7*max(acid)+.3*abs(acid[0]-acid[1])-2.5*(high-acid[0])
                cache.append({'aa':aa,'rotamer':rot,'data':dat,'atoms':atoms,'one':one,'ib':ib,'hs':hshift,'tf':fields,'pre':pre,'center':center})
        # Keep chemical alternatives instead of filling a site with one residue type.
        cache.sort(key=lambda a:a['pre']);chosen=[];counts={}
        for c in cache:
            lim=4 if c['aa'] in 'HDEY' else 2
            if counts.get(c['aa'],0)>=lim:continue
            chosen.append(c);counts[c['aa']]=counts.get(c['aa'],0)+1
            if len(chosen)>=maxoptions:break
        if not chosen:
            # No forced placeholder rescue; reject this backbone explicitly.
            raise ValueError('No geometrically allowed rotamer at loop residue '+str(site+1))
        opts.append(chosen);stats.append({'position':int(site+1),'allowed_rotamers':len(cache),'retained':len(chosen)})
    n=len(opts);m=max(len(a) for a in opts)
    one=np.full((n,m),1e5);ib=np.zeros((n,m,2));hs=ib.copy();tf=np.zeros((n,m,2,len(hpos)));aaidx=np.zeros((n,m),np.int64)
    data=np.zeros((n,m,11,10));sizes=np.zeros((n,m),np.int64);nopt=np.array([len(a) for a in opts])
    for i,options in enumerate(opts):
        for j,c in enumerate(options):
            one[i,j]=c['one'];ib[i,j]=c['ib'];hs[i,j]=c['hs'];tf[i,j]=c['tf'];aaidx[i,j]=AA.index(c['aa']);sizes[i,j]=len(c['data']);data[i,j,:len(c['data'])]=c['data']
    pairs=make_pairs(data,sizes,nopt,bb[sites,1],aaidx)
    return {'bb':bb,'mask':mask,'seq':seq,'orig':orig,'pp':pp,'fixed':fixed,'sites':sites,'opts':opts,'stats':stats,'one':one,'ib':ib,'hs':hs,'tf':tf,'aaidx':aaidx,'pair':pairs,'nopt':nopt,'hpos':hpos}

@njit(cache=True)
def make_pairs(data,sizes,nopt,ca,aaidx):
    n=len(nopt);m=data.shape[1];tab=np.zeros((n,m,n,m),np.float32)
    for i in range(n):
        for j in range(i):
            dist=math.sqrt(np.sum((ca[i]-ca[j])**2))
            if dist>17.:continue
            for u in range(nopt[i]):
                for v in range(nopt[j]):
                    a=data[i,u,:sizes[i,u]];b=data[j,v,:sizes[j,v]]
                    rep,vdw,hb=pair_score(a,b)
                    # Folding regularizer only: no false direct pH binding reward.
                    value=rep+vdw+hb
                    tab[i,u,j,v]=value;tab[j,v,i,u]=value
    return tab

def serialize_design(rec,t,sel,metrics,run,mode):
    seq=list(t['seq']);rows=[];chosen=[]
    bysite={int(s):i for i,s in enumerate(t['sites'])}
    for k in range(len(seq)):
        if k in bysite:
            i=bysite[k];c=t['opts'][i][int(sel[i])];seq[k]=c['aa'];atoms=c['atoms']
            chosen.append({'position':k+1,'aa':c['aa'],'rotamer':int(c['rotamer']),'empirical_source':LIB[c['aa']]['sources'][c['rotamer']]})
        else:atoms=t['fixed'][k]
        rows.append({'position':k+1,'aa':seq[k],'loop':bool(t['mask'][k]),'original_framework_position':int(t['orig'][k]) if not t['mask'][k] else None,'atoms':{a:np.array(p).tolist() for a,p in atoms.items()}})
    # Full source framework ends VTVS; add conventional terminal Ser from same generic fold.
    # Terminal append is deliberately deferred until a coordinate builder exists.
    seq=''.join(seq)
    intervals=[];starts=np.flatnonzero(np.diff(np.r_[False,t['mask'],False]).astype(int))
    for a,b in zip(starts[::2],starts[1::2]):intervals.append([int(a),int(b)])
    cdrs=[seq[a:b] for a,b in intervals]
    return {'candidate_id':rec['pose_id']+'_'+str(run).zfill(2),'pose_id':rec['pose_id'],'sequence':seq,'molecule_class':'nanobody','mode':mode,'seed':20030000+int(rec['pose_id'][1:])*100+run,'metrics':dict(zip(['search_objective','packing_regularizer','human_acid_surrogate','mouse_acid_surrogate','human_neutral_surrogate','ph_energy_contrast_surrogate','histidines_in_loops'],map(float,metrics))),'cdr_intervals_zero_based':intervals,'cdr_sequences':cdrs,'sidechain_sources':chosen,'structure':rows,'status':'unvalidated computational candidate','historical_rng':'MT19937'}

def main(start=0,count=12,nruns=9):
    path=S/'intermediate/docking/poses.json';poses=json.loads(path.read_text());out=S/'intermediate/designs';out.mkdir(exist_ok=True)
    for idx,rec in enumerate(poses[start:start+count],start):
        tick=time.time();pid=rec['pose_id'];log={'pose':pid,'index':idx,'start':time.time()}
        if (out/(pid+'_summary.json')).exists():continue
        try:
            t=candidate_table(rec);results=[];seen=set()
            for run in range(nruns):
                mode=run%3;seed=20030000+int(pid[1:])*100+run
                sel,metrics=anneal(t['nopt'],t['one'],t['pair'],t['ib'],t['hs'],t['tf'],t['aaidx'],seed,6000,mode,6.0 if run>=6 else 6.3)
                d=serialize_design(rec,t,sel,metrics,run,mode)
                if d['sequence'] in seen:continue
                seen.add(d['sequence']);results.append(d)
                (out/(d['candidate_id']+'.json')).write_text(json.dumps(d))
            log.update({'candidate_count':len(results),'site_stats':t['stats'],'seconds':time.time()-tick,'status':'completed','hpos':t['hpos']})
            (out/(pid+'_summary.json')).write_text(json.dumps(log,indent=2))
            print(pid,'candidates',len(results),'sites',len(t['sites']),'seconds',round(time.time()-tick,2),'best',round(min(d['metrics']['search_objective'] for d in results),2),flush=True)
        except Exception as e:
            import traceback;traceback.print_exc();log.update({'status':'failed','error':str(e),'seconds':time.time()-tick});(out/(pid+'_failed.json')).write_text(json.dumps(log,indent=2));print(pid,'FAILED',e,flush=True)
if __name__=='__main__':
    main(*map(int,sys.argv[1:]))
