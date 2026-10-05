"""Classical loop stabilization of our own untested de novo designs.
Finite rotamer enumeration, Gly->Ala/X->Pro entropy rationale (1987), and
geometric disulfide design with Powell minimization (1964). No sequence models.
These are hypotheses; no entropy bonus is equated with measured stability.
"""
from pathlib import Path
import sys,json,copy,itertools,math,hashlib,re
import numpy as np
R=Path('/mnt/data/egfr_campaign');S=R/'stage4';S3=R/'stage3'
sys.path[:0]=[str(S3/'code'),str(R/'stage2/code'),str(R/'code'),str(R/'stage2/runtime/python')]
from classical_design import *
from mm_refine import app,u
from legacy_geometry import dihedral
from scipy.spatial import cKDTree
from scipy.optimize import minimize
from Bio.SeqUtils import seq1


def load(cid):
    src=S3/'intermediate/designs'/f'{cid}.json';parent=json.loads(src.read_text());pfile=S3/'intermediate/refined'/f'{cid}_human6ARU_relaxed.pdb';p=app.PDBFile(str(pfile));xyz=np.asarray(p.positions.value_in_unit(u.angstrom));rows=copy.deepcopy(parent['structure']);rs=[r for r in p.topology.residues() if r.chain.id=='B']
    for row,res in zip(rows,rs):row['atoms']={a.name:xyz[a.index].tolist() for a in res.atoms() if a.element!=app.element.hydrogen and a.name!='OXT'}
    target=[{'aa':seq1(r.name,custom_map={'HIE':'H','HID':'H','HIP':'H'}),'atoms':{a.name:xyz[a.index] for a in r.atoms() if a.element!=app.element.hydrogen}} for r in p.topology.residues() if r.chain.id=='A']
    tenv=np.concatenate([atom_data(r['aa'],r['atoms']) for r in target]);mouse=json.loads((R/'intermediate/mouseAF_aligned.json').read_text());menv=np.concatenate([atom_data(r['aa'],r['atoms']) for r in mouse if 334<=r['human_pos']<=505]);gly=json.loads((R/'stage2/intermediate/target_glycans.json').read_text());gt=cKDTree([r['xyz'] for r in gly]);return parent,rows,tenv,menv,gt,src,pfile


def phi(rows,i):return math.degrees(dihedral(np.array(rows[i-1]['atoms']['C']),*[np.array(rows[i]['atoms'][n]) for n in ['N','CA','C']]))


def rotamer_option(rows,index,aa,tenv,menv,gt):
    bb=np.array([rows[index]['atoms'][n] for n in BB]);env=np.concatenate([atom_data(r['aa'],r['atoms']) for j,r in enumerate(rows) if j!=index]);opts=[]
    ph=phi(rows,index)
    if aa=='P' and not -95<=ph<=-40:return None
    if rows[index]['aa']=='G' and aa!='G' and ph>0:return None
    for rn in range(len(LIB[aa]['xyz'])):
        atoms=conformer(aa,rn,bb);data=atom_data(aa,atoms,LIB[aa]['names']);rp,vw,hb=pair_score(data,env)
        # Proline's N-CD bond must close, not be treated as a favorable clash.
        if aa=='P':
            ring=float(np.linalg.norm(atoms['CD']-atoms['N']))
            if not 1.25<=ring<=1.70:continue
        if rp>20 or (len(data) and gt.query(data[:,:3])[0].min()<2.45):continue
        tv=[]
        for t in [tenv,menv]:
            a,b,c=pair_score(data,t)
            if a>10:break
            tv.append(a+b+c)
        if len(tv)!=2:continue
        val=rp+vw+hb+max(tv);opts.append((val,rn,atoms))
    return min(opts,key=lambda z:z[0]) if opts else None


def rotate_sg(atoms,angle):
    a={k:np.array(v) for k,v in atoms.items()};axis=a['CB']-a['CA'];axis/=np.linalg.norm(axis);v=a['SG']-a['CB'];a['SG']=a['CB']+v*math.cos(angle)+np.cross(axis,v)*math.sin(angle)+axis*np.dot(axis,v)*(1-math.cos(angle));return a


def angle3(a,b,c):
    x=a-b;y=c-b;return math.acos(float(np.clip(np.dot(x,y)/(np.linalg.norm(x)*np.linalg.norm(y)),-1,1)))


def disulfide_option(rows,i,j,tenv,menv,gt):
    if rows[i]['aa']=='G' and phi(rows,i)>0:return None
    if rows[j]['aa']=='G' and phi(rows,j)>0:return None
    env=np.concatenate([atom_data(r['aa'],r['atoms']) for k,r in enumerate(rows) if k not in [i,j]])
    bbs=[np.array([rows[k]['atoms'][n] for n in BB]) for k in [i,j]];best=None
    for ri,rj in itertools.product(range(len(LIB['C']['xyz'])),repeat=2):
        initial=[conformer('C',rn,bb) for rn,bb in zip([ri,rj],bbs)]
        def calculate(x,full=False):
            a,b=[rotate_sg(a,t) for a,t in zip(initial,x)];dis=np.linalg.norm(a['SG']-b['SG']);angles=[angle3(a['CB'],a['SG'],b['SG']),angle3(a['SG'],b['SG'],b['CB'])];chi=float(dihedral(a['CB'],a['SG'],b['SG'],b['CB']))
            geom=100*(dis-2.03)**2+20*sum((v-math.radians(104.))**2 for v in angles)+5*math.cos(chi)**2
            aa=atom_data('C',a,LIB['C']['names']);bb=atom_data('C',b,LIB['C']['names']);val=geom
            for dat in [aa,bb]:
                rp,vw,hb=pair_score(dat,env);val+=rp+vw+hb
                val+=max(sum(pair_score(dat,tenv)),sum(pair_score(dat,menv)))
            if full:return val,a,b,dis,angles,chi,geom
            return val
        opt=minimize(calculate,[0.,0.],method='Powell',options={'maxiter':70,'xtol':1e-4,'ftol':1e-5});res=calculate(opt.x,True);val,a,b,dis,angles,chi,geom=res
        if not (1.9<dis<2.18 and all(abs(math.degrees(v)-104)<25 for v in angles) and abs(abs(math.degrees(chi))-90)<35):continue
        if min(gt.query(np.array([a['CB'],a['SG'],b['CB'],b['SG']]))[0].min() if False else False:pass
        dat=np.array([a['CB'],a['SG'],b['CB'],b['SG']]);mind=float(gt.query(dat)[0].min())
        if mind<2.45:continue
        # No accidental close contact to pre-existing cysteines.
        other=np.array([r['atoms']['SG'] for k,r in enumerate(rows) if r['aa']=='C' and k not in [i,j]])
        if len(other) and np.min(np.linalg.norm(dat[[1,3]][:,None]-other[None,:],axis=2))<3.1:continue
        # Exclude new sulfur's intended partner before scoring other sterics.
        if max(pair_score(atom_data('C',a,LIB['C']['names']),env)[0],pair_score(atom_data('C',b,LIB['C']['names']),env)[0])>25:continue
        if best is None or val<best[0]:best=(val,a,b,{'positions':[i+1,j+1],'rotamers':[ri,rj],'chi_adjustments_radians':opt.x.tolist(),'SG_distance_A':float(dis),'angles_deg':[math.degrees(v) for v in angles],'chi3_deg':math.degrees(chi),'bond_geom_penalty':float(geom),'optimization':'Powell1964'})
    return best


def main():
    out=S/'intermediate/designs';out.mkdir(parents=True,exist_ok=True);records=[];serial=0;seen=set()
    proposals=[{p:'P'} for p in [27,29,35,55,56,57,59,61,64]]+[{p:'A'} for p in [35,51,99,102,111,113]]
    proposals +=[{35:'A',51:'A'},{35:'A',51:'A',111:'A'},{55:'P',51:'A'},{56:'P',51:'A'},{57:'P',51:'A'},{59:'P',51:'A'},{55:'P',59:'P'},{56:'P',64:'P'},{57:'P',64:'P'},{51:'A',99:'A',102:'A',111:'A'},{29:'P',55:'P',111:'A'},{29:'P',57:'P',111:'A'}]
    def save(parent,rows,src,pfile,edits,atoms,score,extra):
        nonlocal serial
        d=copy.deepcopy(parent);seq=list(d['sequence']);newrows=copy.deepcopy(rows);ed=[]
        for position,aa in edits.items():
            i=position-1;old=seq[i];seq[i]=aa;newrows[i]['aa']=aa;newrows[i]['atoms']={n:np.asarray(v).tolist() for n,v in atoms[position].items()};ed.append({'position':position,'from':old,'to':aa})
        seq=''.join(seq)
        if seq in seen:return None
        if re.search(r'N[^P][ST]',seq):return None
        seen.add(seq);serial+=1;cid=f'S{serial:05d}';d.update(candidate_id=cid,sequence=seq,structure=newrows,branch='own_de_novo_loop_stabilization',campaign_stage=4,status='unvalidated computational candidate',seed=None,mutation_selection='deterministic enumeration and Powell optimization')
        d['cdr_sequences']=[seq[a:b] for a,b in d['cdr_intervals_zero_based']];d['ancestry']={'parent_candidate':parent['candidate_id'],'parent_type':'own untested de novo design','parent_file':str(src),'parent_sha256':hashlib.sha256(src.read_bytes()).hexdigest(),'relaxed_backbone_source':str(pfile),'relaxed_backbone_sha256':hashlib.sha256(pfile.read_bytes()).hexdigest()};d['stabilization']={'edits':ed,'local_objective':float(score),'not_a_measured_stability_change':True,**extra};d['metrics']={'scores_valid':False,'packing_regularizer':float(score),'human_acid_surrogate':0.,'mouse_acid_surrogate':0.,'ph_energy_contrast_surrogate':0.,'note':'No affinity or entropy claim; full matched validation required'}
        if 'disulfide' in extra:d['intended_extra_disulfide_pairs']=[extra['disulfide']['positions']]
        d['sidechain_sources']=[s for s in d['sidechain_sources'] if s['position'] not in edits]+extra.get('rotamer_sources',[])
        (out/f'{cid}.json').write_text(json.dumps(d));records.append({'candidate_id':cid,'parent':parent['candidate_id'],'edits':ed,'score':float(score),'extra_disulfide':d.get('intended_extra_disulfide_pairs',[])});print(cid,parent['candidate_id'],ed,'score',round(score,3),flush=True);return cid
    for cid in ['C00003','C00009']:
        parent,rows,te,me,gt,src,pfile=load(cid);cache={}
        for edit in proposals:
            atoms={};score=0.;sources=[];reject=False
            for pos,aa in edit.items():
                key=(pos,aa)
                if key not in cache:cache[key]=rotamer_option(rows,pos-1,aa,te,me,gt)
                opt=cache[key]
                if opt is None:reject=True;break
                score+=opt[0];atoms[pos]=opt[2];sources.append({'position':pos,'aa':aa,'rotamer':opt[1],'source':LIB[aa]['sources'][opt[1]]})
            if reject:records.append({'parent':cid,'proposal':edit,'status':'no compatible rotamer'});continue
            # Reject steric collisions between newly chosen side chains.
            pairmax=0.
            for a,b in itertools.combinations(edit,2):pairmax=max(pairmax,pair_score(atom_data(edit[a],atoms[a],LIB[edit[a]]['names']),atom_data(edit[b],atoms[b],LIB[edit[b]]['names']))[0])
            if pairmax>15:records.append({'parent':cid,'proposal':edit,'status':'new sidechain pair clash'});continue
            save(parent,rows,src,pfile,edit,atoms,score,{'strategy':'Gly-to-Ala and phi-compatible Pro; Matthews1987 rationale only','rotamer_sources':sources})
        for left,right in [(57,60),(36,51),(34,99),(35,73),(108,110)]:
            opt=disulfide_option(rows,left-1,right-1,te,me,gt)
            if opt is None:records.append({'parent':cid,'proposal':{left:'C',right:'C'},'status':'no geometrically acceptable cystine'});continue
            score,a,b,meta=opt;save(parent,rows,src,pfile,{left:'C',right:'C'},{left:a,right:b},score,{'strategy':'geometric additional disulfide hypothesis','disulfide':meta,'rotamer_sources':[{'position':p,'aa':'C','rotamer':rn,'source':LIB['C']['sources'][rn]} for p,rn in zip([left,right],meta['rotamers'])]})
    (S/'reference/stabilization_designs.json').write_text(json.dumps(records,indent=2));print('GENERATED',serial,flush=True)
if __name__=='__main__':main()
