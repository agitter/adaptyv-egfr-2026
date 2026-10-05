"""Finite classical mutations to OUR de novo H358/H383 interface.

Preserve aromatic anchor Y100 and receptor-contacting E103/D54. Test removal
of non-gating positive charges and a second acid contact. No known-binder
sequence is used. Energy decomposition motivates, but does not validate, edits.
"""
from ablation_design import *
from Bio.SeqUtils import seq1


def main():
    cid='B00013';parent_file=S3/'intermediate/designs'/(cid+'.json');parent=json.loads(parent_file.read_text());pdbpath=S3/'intermediate/refined'/(cid+'_human6ARU_relaxed.pdb');pdb=app.PDBFile(str(pdbpath));xyz=np.asarray(pdb.positions.value_in_unit(u.angstrom));rows=copy.deepcopy(parent['structure']);residues=[r for r in pdb.topology.residues() if r.chain.id=='B']
    for row,res in zip(rows,residues):row['atoms']={a.name:xyz[a.index].tolist() for a in res.atoms() if a.element!=app.element.hydrogen and a.name!='OXT'}
    target=[]
    for res in pdb.topology.residues():
        if res.chain.id!='A':continue
        target.append({'aa':seq1(res.name,custom_map={'HIE':'H','HID':'H','HIP':'H'}),'human_pos':int(res.id)+333,'atoms':{a.name:xyz[a.index].tolist() for a in res.atoms() if a.element!=app.element.hydrogen}})
    tenv=np.concatenate([atom_data(r['aa'],r['atoms']) for r in target]);mouse=json.loads((R/'intermediate/mouseAF_aligned.json').read_text());menv=np.concatenate([atom_data(r['aa'],r['atoms']) for r in mouse if 334<=r['human_pos']<=505]);gly=json.loads((R/'stage2/intermediate/target_glycans.json').read_text());tree=cKDTree([r['xyz'] for r in gly])
    edits=[{36:'Q'},{105:'Q'},{36:'Q',105:'Q'},{36:'S'},{105:'S'},{36:'Q',105:'S'},
           {52:'E',36:'Q'},{52:'E',105:'Q'},{52:'E',36:'Q',105:'Q'},
           {53:'D',36:'Q'},{53:'E',36:'Q'},{53:'D',105:'Q'},{53:'E',105:'Q'},
           {52:'E',53:'D'},{52:'E',53:'E'},{52:'D',53:'E'},
           {52:'E',53:'D',36:'Q',105:'Q'},{52:'E',53:'E',36:'Q',105:'Q'},
           {36:'A',105:'Q'},{36:'Q',105:'N'},{36:'Q',105:'Y'},
           {52:'E',36:'Q',105:'N'},{53:'E',36:'Q',105:'Q'},{52:'D',36:'Q',105:'Q'}]
    output=[];n=0
    for edit in edits:
        indices={p-1 for p in edit};fixed=np.concatenate([atom_data(r['aa'],r['atoms']) for i,r in enumerate(rows) if i not in indices]);options={}
        for position,aa in edit.items():
            index=position-1;assert rows[index]['loop'];bb=np.array([rows[index]['atoms'][k] for k in BB]);possible=[]
            for rn in range(len(LIB[aa]['xyz'])):
                atoms=conformer(aa,rn,bb);dat=atom_data(aa,atoms,LIB[aa]['names']);rp,vw,hb=pair_score(dat,fixed)
                if rp>13 or tree.query(dat[:,:3])[0].min()<2.45:continue
                tv=[]
                for t in [tenv,menv]:
                    a,b,c=pair_score(dat,t)
                    if a>10:break
                    tv.append(a+b+c)
                if len(tv)!=2:continue
                score=rp+vw+hb+max(tv);possible.append((float(score),rn,atoms,dat))
            possible.sort(key=lambda r:r[0]);options[position]=possible[:4]
        if any(not v for v in options.values()):output.append({'edit':edit,'status':'no feasible rotamer'});continue
        best=None
        for combo in itertools.product(*options.values()):
            score=sum(o[0] for o in combo)
            for left,right in itertools.combinations(combo,2):
                a,b,c=pair_score(left[3],right[3]);score+=a+b+c
            if best is None or score<best[0]:best=(score,combo)
        n+=1;d=copy.deepcopy(parent);new='C%05d'%n;d['candidate_id']=new;d['structure']=copy.deepcopy(rows);seq=list(d['sequence']);ed=[]
        for (position,aa),opt in zip(edit.items(),best[1]):
            i=position-1;old=seq[i];seq[i]=aa;d['structure'][i]['aa']=aa;d['structure'][i]['atoms']={k:np.asarray(v).tolist() for k,v in opt[2].items()};ed.append({'position':position,'from':old,'to':aa,'rotamer':opt[1],'source':LIB[aa]['sources'][opt[1]]})
        d['sequence']=''.join(seq);d['cdr_sequences']=[d['sequence'][a:b] for a,b in d['cdr_intervals_zero_based']];d['branch']='own_de_novo_parent_charge_balance';d['ancestry']={'parent_candidate':cid,'parent_sha256':hashlib.sha256(parent_file.read_bytes()).hexdigest(),'parent_type':'own untested de novo computational design','relaxed_backbone_source':str(pdbpath),'relaxed_backbone_sha256':hashlib.sha256(pdbpath.read_bytes()).hexdigest()}
        d['charge_balance_edits']=ed;d['sidechain_sources']=[r for r in d['sidechain_sources'] if r['position'] not in edit]+ed;d.pop('focused_edit',None);d['metrics']={'packing_regularizer':float(best[0]),'human_acid_surrogate':0.,'mouse_acid_surrogate':0.,'ph_energy_contrast_surrogate':0.,'scores_valid':False,'note':'Packing-only mutation selection; direct MM and sequence screening still required.'};d['seed']=None;d['mutation_selection']='deterministic finite enumeration';(S3/'intermediate/designs'/(new+'.json')).write_text(json.dumps(d));output.append({'candidate_id':new,'edits':ed,'objective':float(best[0])})
    (S3/'reference/charge_balance_design.json').write_text(json.dumps(output,indent=2))
    for row in output:print(row,flush=True)
if __name__=='__main__':main()
