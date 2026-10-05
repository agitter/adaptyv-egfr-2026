"""Classical focused side-chain redesign around a second receptor histidine.
Only our own newly designed, untested parent is modified. This branch tests
whether H383 engagement can add proton linkage while preserving H358 contacts.
"""
from ablation_design import *
from clamp_geometry import option_features
import clamp_geometry as cg

def main():
    cid='G00117_C04';parent_file=S3/'intermediate/designs'/(cid+'.json');parent=json.loads(parent_file.read_text());pdbpath=S3/'intermediate/refined'/(cid+'_human6ARU_relaxed.pdb');pdb=app.PDBFile(str(pdbpath));xyz=np.asarray(pdb.positions.value_in_unit(u.angstrom));rows=copy.deepcopy(parent['structure']);residues=[r for r in pdb.topology.residues() if r.chain.id=='B'];assert len(rows)==len(residues)
    for row,res in zip(rows,residues):row['atoms']={a.name:xyz[a.index].tolist() for a in res.atoms() if a.element!=app.element.hydrogen and a.name!='OXT'}
    target=[]
    from Bio.SeqUtils import seq1
    for res in pdb.topology.residues():
        if res.chain.id!='A':continue
        aa=seq1(res.name,custom_map={'HIE':'H','HID':'H','HIP':'H'});target.append({'aa':aa,'human_pos':int(res.id)+333,'atoms':{a.name:xyz[a.index].tolist() for a in res.atoms() if a.element!=app.element.hydrogen}})
    cg.HIST=[r for r in target if r['human_pos']==383];ht=cg.HIST[0];center=np.mean([ht['atoms']['ND1'],ht['atoms']['NE2']],0)
    tenv=np.concatenate([atom_data(r['aa'],r['atoms']) for r in target]);mouse=json.loads((R/'intermediate/mouseAF_aligned.json').read_text());menv=np.concatenate([atom_data(r['aa'],r['atoms']) for r in mouse if 334<=r['human_pos']<=505]);gly=json.loads((R/'stage2/intermediate/target_glycans.json').read_text());tree=cKDTree([r['xyz'] for r in gly])
    proposals=[]
    for index in [51,52,53,99]:
        assert rows[index]['loop'];bb=np.array([rows[index]['atoms'][n] for n in BB]);env=np.concatenate([atom_data(r['aa'],r['atoms']) for i,r in enumerate(rows) if i!=index])
        for aa in 'DE':
            opts=[]
            for rn in range(len(LIB[aa]['xyz'])):
                atoms=conformer(aa,rn,bb);data=atom_data(aa,atoms,LIB[aa]['names']);rp,vw,hb=pair_score(data,env)
                if rp>13 or tree.query(data[:,:3])[0].min()<2.45:continue
                tvals=[]
                for t in [tenv,menv]:
                    a,b,c=pair_score(data,t)
                    if a>10:break
                    tvals.append(a+b+c)
                if len(tvals)!=2:continue
                feature=option_features(aa,atoms)[0];cc=charge_center(aa,atoms);distance=float(np.linalg.norm(cc-center))
                val=rp+vw+hb+max(tvals)-3.5*float(feature.sum())-4.*float(min(feature))+screened_charge(-1.,1.,distance,12.)
                opts.append((val,rn,atoms,feature,distance))
            opts.sort(key=lambda o:o[0])
            # Retain contrasting geometries, not just numerically near duplicates.
            chosen=[]
            for o in opts:
                ca=np.array([o[2][n] for n in CHARGE_NAMES[aa]])
                if all(np.linalg.norm(ca-np.array([v[2][n] for n in CHARGE_NAMES[aa]]))>.5 for v in chosen):chosen.append(o)
                if len(chosen)>=3:break
            for opt in chosen:proposals.append((index,aa,opt))
    output=[]
    for n,(index,aa,opt) in enumerate(proposals,1):
        d=copy.deepcopy(parent);new='B%05d'%n;d['candidate_id']=new;d['structure']=copy.deepcopy(rows);old=d['structure'][index]['aa'];d['structure'][index]['aa']=aa;d['structure'][index]['atoms']={k:np.asarray(v).tolist() for k,v in opt[2].items()};seq=list(d['sequence']);seq[index]=aa;d['sequence']=''.join(seq);d['cdr_sequences']=[d['sequence'][a:b] for a,b in d['cdr_intervals_zero_based']];d['branch']='own_de_novo_parent_second_histidine_contact'
        d['ancestry']={'parent_candidate':cid,'parent_sha256':hashlib.sha256(parent_file.read_bytes()).hexdigest(),'parent_type':'own untested de novo computational design','relaxed_backbone_source':str(pdbpath)}
        d['focused_edit']={'position':index+1,'from':old,'to':aa,'rotamer':opt[1],'source':LIB[aa]['sources'][opt[1]],'H383_contact_features':opt[3].tolist(),'distance_to_H383_center_A':opt[4]}
        d['parent_metrics_not_recomputed']=d.pop('metrics');d['metrics']={'packing_regularizer':float(opt[0]),'human_acid_surrogate':0.,'mouse_acid_surrogate':0.,'ph_energy_contrast_surrogate':0.,'scores_valid':False,'note':'Targeted packing objective only; direct all-atom validation required.'}
        d['seed']=None;d['mutation_selection']='deterministic finite enumeration';d['sidechain_sources']=[r for r in d['sidechain_sources'] if r['position']!=index+1]+[d['focused_edit']];d.pop('clamp_geometry',None);(S3/'intermediate/designs'/(new+'.json')).write_text(json.dumps(d));output.append({'candidate_id':new,**d['focused_edit'],'local_search_objective':float(opt[0])})
    (S3/'reference/focused_mutations.json').write_text(json.dumps(output,indent=2))
    for o in output:print(o['candidate_id'],o['position'],o['from'],o['to'],o['rotamer'],round(o['distance_to_H383_center_A'],2),[round(x,3) for x in o['H383_contact_features']],round(o['local_search_objective'],2))
if __name__=='__main__':main()
