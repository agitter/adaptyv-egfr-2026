"""Verify every new candidate was searched in every reference database."""
from pathlib import Path
import json,hashlib
R=Path('/mnt/data/egfr_campaign');S=R/'stage5';O=S/'intermediate/blast'

def main():
 catalog=json.loads((S/'reference/candidate_catalog.json').read_text());expected={r['candidate_id'] for r in catalog};cdr=[]
 for name in ['stage5_cdr3_screen.json','tail_cdr3_screen.json','chemistry_cdr3_screen.json']:cdr+=json.loads((S/'intermediate/novelty'/name).read_text())
 assert len(cdr)==len(expected) and {r['candidate_id'] for r in cdr}==expected
 summaries=[];per={cid:{} for cid in expected}
 for db in ['swissprot','pdb','antibodies_augmented']:
  observed=[];hits={};sources=[];controls=[]
  for group,tag in [('stage5',''),('tail','_tail'),('chemistry','_chemistry')]:
   summary=S/'reference'/f'blast_{db}{tag}_summary.json';d=json.loads(summary.read_text());p=O/f'{db}_{group}.tsv';query=[l.split(':',1)[1].strip().split()[0] for l in p.read_text().splitlines() if l.startswith('# Query:')];ids=[q for q in query if q in expected];observed+=ids;hits.update({k:v for k,v in d['hits'].items() if k in expected});assert d['all_query_ids_verified'] and d['positive_control_exact_match'];assert len(ids)==d['query_candidates'];sources.append({'summary':str(summary),'summary_sha256':hashlib.sha256(summary.read_bytes()).hexdigest(),'results':str(p),'results_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'candidate_queries':len(ids)});controls.append({'group':group,'exact_control_pass':True,'source':d.get('positive_control_source','see resumable-search original query/provenance')})
  assert len(observed)==len(expected) and set(observed)==expected and len(set(observed))==len(observed)
  maximum=[]
  for cid in expected:
   cover=[r for r in hits.get(cid,[]) if r['query_coverage_pct']>=90];val=max([r['identity_pct'] for r in cover],default=None);per[cid][db]={'search_complete':True,'reported_hit_count':len(hits.get(cid,[])),'maximum_reported_identity_at_90pct_query_coverage':val}
   if val is not None:maximum.append(val)
  summaries.append({'database':db,'candidate_queries':len(observed),'all_ids_verified':True,'maximum_reported_identity_pct_at_90pct_coverage':max(maximum,default=None),'positive_controls':controls,'sources':sources})
 cdrmap={r['candidate_id']:r for r in cdr};rows=[]
 for item in catalog:
  cid=item['candidate_id'];r=cdrmap[cid];a=json.loads((S/'intermediate/evaluation'/f'{cid}.json').read_text());rows.append({'candidate_id':cid,'local_initial_gate_pass':a['pass'],'cdr3_edit_identity':r['max_edit_identity'],'cdr3_local_novelty_pass':r['below_70pct'],'initial_and_local_novelty_pass':a['pass'] and r['below_70pct'],'whole_chain_searches':per[cid]})
 result={'new_candidates':len(rows),'unique_CDR3s':len({r['cdr3_query'] for r in cdr}),'CDR3_reference_count':203818,'cdr3_local_pass':sum(r['below_70pct'] for r in cdr),'cdr3_rejected_ids':[r['candidate_id'] for r in cdr if not r['below_70pct']],'maximum_CDR3_edit_identity':max(r['max_edit_identity'] for r in cdr),'initial_gate_pass':sum(r['local_initial_gate_pass'] for r in rows),'initial_and_local_novelty_pass':sum(r['initial_and_local_novelty_pass'] for r in rows),'databases':summaries,'records':rows,'qualification':'All three local searches complete; local CDR3 motif/normalization is not official IMGT annotation. BLAST top50 hits and these database snapshots are not the organizers full novelty pipeline or upload acceptance. Shared antibody framework similarity is expected.'};(S/'reference/novelty_complete.json').write_text(json.dumps(result,indent=2));print({k:v for k,v in result.items() if k not in ['records','databases']});print([(d['database'],d['maximum_reported_identity_pct_at_90pct_coverage']) for d in summaries])
if __name__=='__main__':main()
