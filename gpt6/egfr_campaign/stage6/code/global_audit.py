"""Uniform classical geometry and independent-prior audit of saved models.
These calculations do not establish folding, affinity or neutral-pH nondetection.
"""
from pathlib import Path
import sys,json,math,hashlib,itertools,time,re
import numpy as np
from scipy.special import logsumexp
from scipy.spatial import cKDTree
R=Path('/mnt/data/egfr_campaign');S=R/'stage6'
sys.path[:0]=[str(R/x/'code') for x in ['stage5','stage4','stage3','stage2']]+[str(R/'code'),str(R/'stage2/runtime/python')]
from evaluate_pool import structural_audit,sequence_audit,TARGETS
from audit_refined import audit as strict_audit
from mm_refine import app,u
from proton_ensemble import RT
AA=set('ACDEFGHIKLMNPQRSTVWY')
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def independent_priors(d):
 micro=d['microstates'];states=np.asarray([x['state'] for x in micro],int);n=states.shape[1]
 if n>4:return {'status':'not expanded: more than four sites'}
 energies=np.asarray([[p['delta_kcal_proxy'] for p in m['energy']['states']] for m in micro]);opts=np.array(list(itertools.product([5.8,6.3,6.8],[.2,.5,.8])));pr=np.array(list(itertools.product(range(9),repeat=n)));p=opts[pr,0];f=opts[pr,1];out=[]
 for ph in [6.5,7.4]:
  lw=math.log(10)*(p-ph)@(states==2).T+np.log(f)@(states==0).T+np.log(1-f)@(states==1).T;lzf=logsumexp(lw,axis=1)
  out.append(np.stack([-RT*(logsumexp(lw-energies[:,j][None,:]/RT,axis=1)-lzf) for j in range(energies.shape[1])],axis=1))
 contrast=out[1]-out[0];assert np.max(np.abs(contrast))<=n*.9*RT*math.log(10)+1e-8;pi,mi=np.unravel_index(contrast.argmin(),contrast.shape);err=[]
 for old in d['proton_model']['scenarios']:
  ix=np.flatnonzero((p==old['assumed_unbound_pKa']).all(1)&(f==old['assumed_unbound_HIE_fraction']).all(1));js=[j for j,v in enumerate(micro[0]['energy']['states']) if v['solute_dielectric']==old['solute_dielectric'] and v['kappa_nm_inverse']==old['kappa_nm_inverse']];assert len(ix)==len(js)==1;err.append(abs(contrast[ix[0],js[0]]-old['contrast_kcal']))
 assert max(err)<1e-9
 return {'status':'complete','minimum_contrast_kcal':float(contrast.min()),'maximum_contrast_kcal':float(contrast.max()),'scenario_count':int(contrast.size),'all_acid_on':bool(np.all(contrast>0)),'tied_subset_max_error_kcal':float(max(err)),'worst_case':{'pKas':p[pi].tolist(),'HIE_fractions':f[pi].tolist(),'physical':micro[0]['energy']['states'][mi],'acid_proxy':float(out[0][pi,mi]),'neutral_proxy':float(out[1][pi,mi])},'qualification':'Assumed independent priors; fixed-coordinate proxy, not calibrated affinity or nondetection.'}
def model_audit(record):
 dp=Path(record['design_source']);d=json.loads(dp.read_text());path=Path(record['coordinate_source']);p=app.PDBFile(str(path));xyz=np.asarray(p.positions.value_in_unit(u.angstrom));binder=[];target={}
 for res in p.topology.residues():
  atoms={a.name:xyz[a.index].tolist() for a in res.atoms() if a.element!=app.element.hydrogen and a.name!='OXT'}
  if res.chain.id=='B':binder.append(atoms)
  elif res.chain.id=='A':target[int(res.id)+333]=atoms
 assert len(binder)==len(d['structure'])
 for row,atoms in zip(d['structure'],binder):row['atoms']=atoms
 base=structural_audit(d);species='mouseAF' if record['species']=='mouse' else 'human6ARU';whole=np.array([v for row in TARGETS[species] for v in target.get(row['human_pos'],row['atoms']).values()]);bx=np.array([v for row in binder for v in row.values()]);dist=cKDTree(whole).query(bx)[0];strict=strict_audit(path,d)
 ap=path.with_name(path.name.replace('_relaxed.pdb','_audit.json'));old=json.loads(ap.read_text()) if ap.exists() else None
 gate=(strict['pass'] and int(sum(dist<2))==0 and base['resolved_glycans']['atoms_below_2A']==0 and all(base[c]['heavy_atoms_below_2A']==0 for c in ['1IVO','1NQL']) and (old is None or old.get('pass',False)))
 return {**record,'source_sha256':sha(path),'design_sha256':sha(dp),'strict_geometry':strict,'reference_geometry':base,'same_species_relaxed_context':{'minimum_A':float(dist.min()),'atoms_below2A':int(sum(dist<2))},'sequence_audit':sequence_audit(d),'minimizer_audit_present':old is not None,'minimizer_geometry_pass':old.get('pass') if old else None,'combined_geometry_pass':bool(gate),'independent_priors':independent_priors(json.loads(Path(record['model']).read_text()))}
def main():
 start=time.time();models=json.loads((S/'reference/model_inventory.json').read_text());out=S/'intermediate/uniform_models';out.mkdir(exist_ok=True);summ=[]
 for i,row in enumerate(models):
  path=out/(row['stage']+'_'+Path(row['model']).stem+'.json')
  if path.exists():result=json.loads(path.read_text())
  else:result=model_audit(row);path.write_text(json.dumps(result,indent=2))
  small={k:v for k,v in result.items() if k not in ['strict_geometry','reference_geometry','sequence']};small['result_source']=str(path);summ.append(small)
  if i%10==0:print('Audited model',i+1,len(models),round(time.time()-start,1),flush=True)
 (S/'reference/uniform_model_summary.json').write_text(json.dumps(summ,indent=2));novel={}
 for stage,names in [('stage2',['candidate_cdr3_screen_augmented.json']),('stage3',['candidate_cdr3_screen_augmented.json']),('stage4',['stage4_cdr3_screen.json']),('stage5',['stage5_cdr3_screen.json','tail_cdr3_screen.json','chemistry_cdr3_screen.json'])]:
  for name in names:
   p=R/stage/'intermediate/novelty'/name
   for x in json.loads(p.read_text()):novel[stage+':'+x['candidate_id']]={**x,'source':str(p)}
 bycdr={}
 for n in novel.values():
  c=n.get('cdr3_query');old=bycdr.get(c)
  if old is None or n['max_edit_identity']>old['max_edit_identity']:bycdr[c]=n
 allrows=[];byseq={}
 for m in summ:byseq.setdefault(m['sequence_sha256'],[]).append(m)
 for stage in ['stage2','stage3','stage4','stage5']:
  for path in sorted((R/stage/'intermediate/designs').glob('*.json')):
   d=json.loads(path.read_text())
   if not isinstance(d,dict) or 'sequence' not in d:continue
   key=stage+':'+d['candidate_id'];seq=d['sequence'];assert seq==''.join(r['aa'] for r in d['structure']) and len(seq)==len(d['structure']) and set(seq)<=AA and 10<=len(seq)<=250,key
   g=structural_audit(d);q=sequence_audit(d);nv=novel.get(key,bycdr.get(d['cdr_sequences'][2]));mods=byseq.get(hashlib.sha256(seq.encode()).hexdigest(),[]);initial=g['geometry_prefilter_pass'] and all(g[c]['heavy_atoms_below_2A']==0 for c in ['1IVO','1NQL']) and not q['liabilities']['glycosylation_sequons']
   row={'key':key,'sequence':seq,'design_source':str(path),'design_sha256':sha(path),'length':len(seq),'branch':d.get('branch','initial'),'pose_id':d.get('pose_id'),'CDRs':d['cdr_sequences'],'intervals':d['cdr_intervals_zero_based'],'sequence_audit':q,'initial_geometry':g,'initial_context_pass':bool(initial),'local_novelty_pass':bool(nv is not None and nv['max_edit_identity']<.7),'local_novelty':nv,'model_keys':[m['model'] for m in mods],'any_relaxed_human_pass':any(m['species']=='human' and m['combined_geometry_pass'] for m in mods),'any_relaxed_mouse_pass':any(m['species']=='mouse' and m['combined_geometry_pass'] for m in mods),'initial_metrics':d.get('metrics',{})}
   ev=R/stage/'intermediate/evaluation'/path.name
   if ev.exists():
    e=json.loads(ev.read_text());row['coarse_evaluation_source']=str(ev);row['coarse_human_proton_model']=e.get('human_proton_model');row['coarse_mouse_proton_model']=e.get('mouse_proton_model')
   allrows.append(row)
 (S/'reference/campaign_audit.json').write_text(json.dumps(allrows,indent=2));counts={'records':len(allrows),'unique_sequences':len(set(x['sequence'] for x in allrows)),'uniform_models':len(summ),'uniform_model_geometry_pass':sum(m['combined_geometry_pass'] for m in summ),'initial_context_and_local_novelty':sum(x['initial_context_pass'] and x['local_novelty_pass'] for x in allrows),'seconds':time.time()-start};(S/'reference/campaign_audit_summary.json').write_text(json.dumps(counts,indent=2));print(counts,flush=True)
if __name__=='__main__':main()
