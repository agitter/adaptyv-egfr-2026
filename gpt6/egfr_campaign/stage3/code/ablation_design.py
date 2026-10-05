"""Targeted redesign of OUR untested de novo stage-2 candidate, not a known binder.

The full three-state diagnostic identified adverse binder-His protonation at
positions 108/109. This branch neutralizes these sites while retaining the
independently created receptor-His contacts of the parent. Classical rotamer
packing is repeated, with the parent side-chain coordinates otherwise fixed.
"""
from pathlib import Path
import json,copy,sys,itertools,hashlib
import numpy as np
R=Path('/mnt/data/egfr_campaign');S3=R/'stage3';sys.path[:0]=[str(R/'stage2/code'),str(R/'code'),str(R/'stage2/runtime/python')]
from classical_design import *
from mm_refine import app,u
from scipy.spatial import cKDTree


def main():
    parent_file=R/'stage2/intermediate/designs/P00024_07.json';parent=json.loads(parent_file.read_text());pdbpath=R/'stage2/intermediate/mm_binding/P00024_07_human6ARU_complex_relaxed.pdb';pdb=app.PDBFile(str(pdbpath));xyz=np.asarray(pdb.positions.value_in_unit(u.angstrom));rows=copy.deepcopy(parent['structure']);residues=[r for r in pdb.topology.residues() if r.chain.id=='B'];assert len(rows)==len(residues)
    for row,res in zip(rows,residues):row['atoms']={a.name:xyz[a.index].tolist() for a in res.atoms() if a.element!=app.element.hydrogen and a.name!='OXT'}
    pos=[107,108];assert all(rows[i]['aa']=='H' for i in pos)
    fixed=np.concatenate([atom_data(row['aa'],row['atoms']) for i,row in enumerate(rows) if i not in pos])
    targets=[]
    for species in ['human6ARU','mouseAF']:
        rs=json.loads((R/'intermediate'/(species+'_aligned.json')).read_text());targets.append(np.concatenate([atom_data(r['aa'],r['atoms']) for r in rs if 334<=r['human_pos']<=505]))
    gly=json.loads((R/'stage2/intermediate/target_glycans.json').read_text());tree=cKDTree([r['xyz'] for r in gly])
    allowed='NQAYFST';opts={}
    for i in pos:
        bb=np.array([rows[i]['atoms'][n] for n in BB]);options={}
        for aa in allowed:
            possible=[]
            for rn in range(len(LIB[aa]['xyz'])):
                atoms=conformer(aa,rn,bb);data=atom_data(aa,atoms,LIB[aa]['names']);clash=tree.query(data[:,:3])[0].min() if len(data) else 99.
                if clash<2.45:continue
                rp,vw,hb=pair_score(data,fixed)
                if rp>15:continue
                val=rp+vw+hb;target=[]
                for t in targets:
                    r,v,h=pair_score(data,t)
                    if r>12:break
                    target.append(r+v+h)
                if len(target)!=2:continue
                val+=max(target);possible.append((float(val),rn,atoms,data))
            possible.sort(key=lambda r:r[0]);options[aa]=possible[:4]
        opts[i]=options
    output=[];n=0
    pairs=['QQ','NN','QN','NQ','YY','YQ','QY','FQ','QA','AA','TS','SY']
    for a,b in pairs:
        best=None
        for left,right in itertools.product(opts[pos[0]][a],opts[pos[1]][b]):
            rep,vw,hb=pair_score(left[3],right[3]);score=left[0]+right[0]+rep+vw+hb
            if best is None or score<best[0]:best=(score,left,right)
        if best is None:output.append({'substitutions':a+b,'status':'no feasible rotamer'});continue
        n+=1;d=copy.deepcopy(parent);cid='A%05d'%n;d['candidate_id']=cid;d['structure']=copy.deepcopy(rows);seq=list(d['sequence']);edits=[]
        for i,aa,opt in zip(pos,[a,b],best[1:]):
            seq[i]=aa;d['structure'][i]['aa']=aa;d['structure'][i]['atoms']={k:np.asarray(v).tolist() for k,v in opt[2].items()};edits.append({'position':i+1,'from':'H','to':aa,'rotamer':opt[1],'source':LIB[aa]['sources'][opt[1]]})
        d['sequence']=''.join(seq);d['cdr_sequences']=[d['sequence'][a:b] for a,b in d['cdr_intervals_zero_based']];d['branch']='own_de_novo_parent_histidine_ablation';d['campaign_stage']=3
        d['ancestry']={'parent_candidate':'P00024_07','parent_source':str(parent_file),'parent_sha256':hashlib.sha256(parent_file.read_bytes()).hexdigest(),'relaxed_backbone_source':str(pdbpath),'criterion':'remove protonation states disfavored by direct MM diagnostic; parent is our own untested computational sequence, not an existing binder'}
        d['sidechain_sources']=[r for r in d['sidechain_sources'] if r['position'] not in [i+1 for i in pos]]+edits;d['ablation_edits']=edits;d['parent_metrics_not_recomputed']=d.pop('metrics')
        # Placeholder only for compatibility with the geometry-evaluation script;
        # these zeros are explicitly not new search-affinity scores.
        d['metrics']={'packing_regularizer':float(best[0]),'human_acid_surrogate':0.,'mouse_acid_surrogate':0.,'ph_energy_contrast_surrogate':0.,'note':'Parent search score invalid after mutation; all-atom validation required.'}
        d['seed']=None;d['mutation_selection']='deterministic finite enumeration';d['status']='unvalidated computational redesign';(S3/'intermediate/designs'/(cid+'.json')).write_text(json.dumps(d));output.append({'candidate_id':cid,'substitutions':a+b,'status':'generated','packing_objective':float(best[0])})
    (S3/'reference/ablation_design.json').write_text(json.dumps(output,indent=2));print(output)
if __name__=='__main__':main()
