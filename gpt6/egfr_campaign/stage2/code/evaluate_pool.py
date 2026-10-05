"""Independent screening of completed sequence/structure candidates."""
from pathlib import Path
import json,time,re,sys,math
import numpy as np
from scipy.spatial import cKDTree
from classical_design import *
from proton_polynomial import evaluate_design

TARGETS={n:json.loads((O/(n+'_aligned.json')).read_text()) for n in ['human6ARU','mouseAF','1NQL','1IVO']}
TREES={n:cKDTree(np.array([p for r in a for p in r['atoms'].values()])) for n,a in TARGETS.items()}
GLY=json.loads((O/'human6ARU_glycans.json').read_text());GTREE=cKDTree(np.array([a['xyz'] for a in GLY]))


def sequence_audit(d):
    seq=d['sequence'];regions=d['cdr_intervals_zero_based'];cdr=''.join(seq[a:b] for a,b in regions)
    liabilities={'glycosylation_sequons':[m.start()+1 for m in re.finditer(r'N[^P][ST]',seq)],'NG_motifs':[m.start()+1 for m in re.finditer(r'NG',seq)],'DG_motifs':[m.start()+1 for m in re.finditer(r'DG',seq)],'unpaired_cdr_cysteines':cdr.count('C'),'max_repeat':max(len(m.group(0)) for m in re.finditer(r'(.)\1*',seq)),'hydrophobic_runs':[m.group(0) for m in re.finditer(r'[VILMFWY]{4,}',seq)],'net_formal_charge_without_histidine':seq.count('K')+seq.count('R')-seq.count('D')-seq.count('E')}
    hard=(10<=len(seq)<=250 and set(seq)<=set(AA) and len(regions)==3 and seq.count('C')==2 and cdr.count('C')==0)
    return {'length':len(seq),'format_pass':hard,'liabilities':liabilities,'cdr_lengths':[b-a for a,b in regions]}


def structural_audit(d):
    rows=d['structure'];xyz=np.array([p for r in rows for p in r['atoms'].values()]);loopxyz=np.array([p for r in rows if r['loop'] for p in r['atoms'].values()]);bb=np.array([p for r in rows for n,p in r['atoms'].items() if n in BB])
    out={}
    for name,tree in TREES.items():
        dis,_=tree.query(xyz);db,_=tree.query(bb)
        out[name]={'minimum_heavy_distance_A':float(dis.min()),'heavy_atoms_below_2A':int((dis<2.).sum()),'minimum_backbone_distance_A':float(db.min())}
    dg,_=GTREE.query(xyz);out['resolved_glycans']={'min_heavy_distance_A':float(dg.min()),'atoms_below_2A':int((dg<2.).sum())}
    # Native stereochemistry is checked without reference to optimizer score.
    chir=[];backbone=[]
    for r in rows:
        a=r['atoms']
        if 'CB' in a:
            v=np.linalg.det(np.array([np.array(a[x])-a['CA'] for x in ['N','C','CB']]))
            if v<=0:chir.append(r['position'])
    for a,b in zip(rows,rows[1:]):
        dist=float(np.linalg.norm(np.array(a['atoms']['C'])-b['atoms']['N']))
        if not 1.15<dist<1.6:backbone.append([a['position'],dist])
    out['alpha_inversions']=chir;out['peptide_bond_outliers']=backbone
    out['geometry_prefilter_pass']=not chir and not backbone and all(out[n]['heavy_atoms_below_2A']==0 for n in ['human6ARU','mouseAF']) and out['resolved_glycans']['atoms_below_2A']==0
    return out


def contact_audit(d):
    contacts=[];rows=d['structure']
    for b in rows:
        if b['aa'] not in 'HDE':continue
        bn=CHARGE_NAMES[b['aa']]
        if not all(k in b['atoms'] for k in bn):continue
        for t in TARGETS['human6ARU']:
            if not ((b['aa']=='H' and t['aa'] in 'DE') or (b['aa'] in 'DE' and t['aa']=='H')):continue
            tn=CHARGE_NAMES[t['aa']]
            if not all(k in t['atoms'] for k in tn):continue
            dis=np.linalg.norm(np.array([b['atoms'][k] for k in bn])[:,None,:]-np.array([t['atoms'][k] for k in tn])[None,:,:],axis=2)
            if dis.min()<4.5:contacts.append({'binder_position':b['position'],'binder_aa':b['aa'],'target_human_position':t['human_pos'],'target_aa':t['aa'],'min_NO_distance_A':float(dis.min()),'target_conserved':bool(t['conserved'])})
    return contacts


def main(start=0,count=10000):
    inp=S/'intermediate/designs';out=S/'intermediate/evaluation';out.mkdir(exist_ok=True)
    files=sorted(p for p in inp.glob('*.json') if '_summary' not in p.name and '_failed' not in p.name)
    t0=time.time()
    for i,p in enumerate(files[start:start+count],start):
        dest=out/p.name
        if dest.exists():continue
        d=json.loads(p.read_text());a=sequence_audit(d);g=structural_audit(d);c=contact_audit(d)
        eval={'candidate_id':d['candidate_id'],'sequence_audit':a,'structure_audit':g,'candidate_ion_pairs':c,'direct_conserved_titration_contacts':len({(x['binder_position'],x['target_human_position']) for x in c if x['target_conserved']})}
        if a['format_pass'] and g['geometry_prefilter_pass']:
            try:
                eval['human_proton_model']=evaluate_design(d['structure'],TARGETS['human6ARU'],True)
                eval['mouse_proton_model']=evaluate_design(d['structure'],TARGETS['mouseAF'],False)
            except Exception as e:eval['proton_model_error']=repr(e)
        dest.write_text(json.dumps(eval))
        if i%25==0:print('evaluated',i+1,'seconds',round(time.time()-t0,1),flush=True)
if __name__=='__main__':main(*map(int,sys.argv[1:]))
