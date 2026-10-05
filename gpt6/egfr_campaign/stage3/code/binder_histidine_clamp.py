"""Classical search for a third proton-linked contact in OUR C00003 design.
A binder histidine must present BOTH ring nitrogens toward target carboxylate
acceptors. Merely placing histidine near any acidic atom is not sufficient.
All side chains are finite empirical rotamers; no learned design model is used.
"""
from ablation_design import *
from Bio.SeqUtils import seq1


def main():
    cid='C00003';source=S3/'intermediate/designs'/(cid+'.json');parent=json.loads(source.read_text());pdbpath=S3/'intermediate/refined'/(cid+'_human6ARU_relaxed.pdb');pdb=app.PDBFile(str(pdbpath));xyz=np.asarray(pdb.positions.value_in_unit(u.angstrom));rows=copy.deepcopy(parent['structure']);br=[r for r in pdb.topology.residues() if r.chain.id=='B']
    for row,res in zip(rows,br):row['atoms']={a.name:xyz[a.index].tolist() for a in res.atoms() if a.element!=app.element.hydrogen and a.name!='OXT'}
    target=[]
    for res in pdb.topology.residues():
        if res.chain.id!='A':continue
        target.append({'aa':seq1(res.name,custom_map={'HIE':'H','HID':'H','HIP':'H'}),'human_pos':int(res.id)+333,'atoms':{a.name:xyz[a.index].tolist() for a in res.atoms() if a.element!=app.element.hydrogen}})
    acceptors=[(r['human_pos'],name,np.asarray(r['atoms'][name])) for r in target if r['aa'] in 'DE' for name in (['OD1','OD2'] if r['aa']=='D' else ['OE1','OE2'])];acxyz=np.array([r[2] for r in acceptors]);tenv=np.concatenate([atom_data(r['aa'],r['atoms']) for r in target]);mouse=json.loads((R/'intermediate/mouseAF_aligned.json').read_text());menv=np.concatenate([atom_data(r['aa'],r['atoms']) for r in mouse if 334<=r['human_pos']<=505]);gly=json.loads((R/'stage2/intermediate/target_glycans.json').read_text());tree=cKDTree([r['xyz'] for r in gly]);options=[]
    for i,row in enumerate(rows):
        if not row['loop'] or i+1 in [54,100,103]:continue
        bb=np.array([row['atoms'][n] for n in BB]);fixed=np.concatenate([atom_data(r['aa'],r['atoms']) for j,r in enumerate(rows) if j!=i])
        for rn in range(len(LIB['H']['xyz'])):
            atoms=conformer('H',rn,bb);data=atom_data('H',atoms,LIB['H']['names']);rp,vw,hb=pair_score(data,fixed)
            if rp>13 or tree.query(data[:,:3])[0].min()<2.45:continue
            tv=[]
            for t in [tenv,menv]:
                a,b,c=pair_score(data,t)
                if a>10:break
                tv.append(a+b+c)
            if len(tv)!=2:continue
            scores=[];contacts=[]
            for name,left,right in [('ND1','CG','CE1'),('NE2','CE1','CD2')]:
                n=np.asarray(atoms[name]);v1=np.asarray(atoms[left])-n;v2=np.asarray(atoms[right])-n;direction=-(v1/np.linalg.norm(v1)+v2/np.linalg.norm(v2));direction/=np.linalg.norm(direction);vec=acxyz-n;ds=np.linalg.norm(vec,axis=1);cs=(vec@direction)/ds;feature=np.exp(-((ds-2.85)/.55)**2)*np.maximum(0.,cs)**2;feature[ds<2.3]=0.;j=int(np.argmax(feature));scores.append(float(feature[j]));contacts.append({'nitrogen':name,'target_position':acceptors[j][0],'target_atom':acceptors[j][1],'distance_A':float(ds[j]),'direction_cosine':float(cs[j])})
            minfeat=min(scores);total=rp+vw+hb+max(tv)-10*minfeat
            if minfeat>.005:options.append((minfeat,total,i,rn,atoms,scores,contacts))
    options.sort(key=lambda r:(-r[0],r[1]));chosen=[];perpos={}
    for o in options:
        if perpos.get(o[2],0)>=2:continue
        chosen.append(o);perpos[o[2]]=perpos.get(o[2],0)+1
        if len(chosen)>=10:break
    output=[]
    for n,o in enumerate(chosen,1):
        feature,score,i,rn,atoms,feats,contacts=o;d=copy.deepcopy(parent);new='D%05d'%n;d['candidate_id']=new;d['structure']=copy.deepcopy(rows);old=rows[i]['aa'];d['structure'][i]['aa']='H';d['structure'][i]['atoms']={k:np.asarray(v).tolist() for k,v in atoms.items()};seq=list(d['sequence']);seq[i]='H';d['sequence']=''.join(seq);d['cdr_sequences']=[d['sequence'][a:b] for a,b in d['cdr_intervals_zero_based']];d['branch']='own_de_novo_parent_binder_His_clamp';d['ancestry']={'parent_candidate':cid,'parent_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'parent_type':'own untested de novo design','relaxed_backbone_source':str(pdbpath),'relaxed_backbone_sha256':hashlib.sha256(pdbpath.read_bytes()).hexdigest()};edit={'position':i+1,'from':old,'to':'H','rotamer':rn,'source':LIB['H']['sources'][rn]};d['histidine_clamp_edit']=edit;d['histidine_clamp_features']=feats;d['histidine_clamp_contacts']=contacts;d['sidechain_sources']=[r for r in d['sidechain_sources'] if r['position']!=i+1]+[edit];d['metrics']={'packing_regularizer':float(score),'human_acid_surrogate':0.,'mouse_acid_surrogate':0.,'ph_energy_contrast_surrogate':0.,'scores_valid':False,'note':'Third-contact geometry proposal, not binding validation.'};d['seed']=None;d['mutation_selection']='deterministic finite enumeration';(S3/'intermediate/designs'/(new+'.json')).write_text(json.dumps(d));output.append({'candidate_id':new,'min_two_nitrogen_feature':feature,'edit':edit,'contacts':contacts})
    record={'parent':cid,'options_with_two_features_above_0p005':len(options),'selected':output,'note':'Weak geometric hypotheses are retained as rejection controls; all require explicit proton ensemble evaluation.'};(S3/'reference/binder_histidine_clamp.json').write_text(json.dumps(record,indent=2));print(json.dumps(record,indent=2),flush=True)
if __name__=='__main__':main()
