"""Verify complete BLAST query coverage and preserve full-chain hit summaries.
Full-chain identities are NOT an antibody novelty verdict; see CDR-specific data.
"""
from pathlib import Path
import json,hashlib
R=Path('/mnt/data/egfr_campaign');S=R/'stage4';D=S/'intermediate/blast'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def parse(path):
 queries=set();hits={}
 for line in path.read_text().splitlines():
  if line.startswith('# Query: '):queries.add(line.split(': ',1)[1].split()[0])
  elif line and not line.startswith('#'):
   z=line.split('\t');assert len(z)==13,(path,line)
   r={'subject':z[1],'percent_identity':float(z[2]),'alignment_length':int(z[3]),'query_length':int(z[4]),'subject_length':int(z[5]),'percent_query_coverage_hsp':int(z[6]),'evalue':float(z[7]),'bitscore':float(z[8]),'query_start':int(z[9]),'query_end':int(z[10]),'subject_start':int(z[11]),'subject_end':int(z[12])}
   hits.setdefault(z[0],[]).append(r)
 return queries,hits

def main():
 queries={p.stem for p in (S/'intermediate/designs').glob('*.json')};assert len(queries)==63
 main=json.loads((S/'reference/blast_stage4_status.json').read_text());follow=json.loads((S/'reference/blast_followup_status.json').read_text());assert main['complete'] and all(r['positive_control_exact_full_match_pass'] for r in follow)
 out={'candidate_count':len(queries),'databases':{},'scope':'BLAST is a heuristic top-20 hit search, not exhaustive global alignment. Generic framework similarities are expected for nanobodies. CDR screening is separate; official classification/annotation and patent coverage are not replicated.'}
 for name in ['swissprot','pdb','antibodies_augmented']:
  a=D/(name+'_stage4.tsv');b=D/(name+'_followup.tsv');q1,h1=parse(a);q2,h2=parse(b);observed=(q1|q2)&queries;assert observed==queries,(name,queries-observed)
  hh={**h1,**h2};records={}
  for cid in sorted(queries):
   h=hh.get(cid,[]);assert all(x['query_length']==125 for x in h)
   broad=[x for x in h if x['percent_query_coverage_hsp']>=90]
   records[cid]={'hit_count':len(h),'highest_bitscore_hit':max(h,key=lambda r:r['bitscore']) if h else None,'highest_identity_hit_coverage90plus':max(broad,key=lambda r:r['percent_identity']) if broad else None}
  maxima=[r['highest_identity_hit_coverage90plus']['percent_identity'] for r in records.values() if r['highest_identity_hit_coverage90plus']]
  out['databases'][name]={'query_count_verified':len(observed),'input_result_hashes':{a.name:sha(a),b.name:sha(b)},'candidate_queries_without_reported_hit':sum(not r['hit_count'] for r in records.values()),'maximum_reported_identity_among_hits_with90plus_query_coverage':max(maxima) if maxima else None,'records':records,'positive_control_pass':True}
 (S/'intermediate/novelty/blast_coverage_and_hits.json').write_text(json.dumps(out,indent=2))
 summary={'candidate_count':63,'database_query_counts':{k:v['query_count_verified'] for k,v in out['databases'].items()},'full_chain_positive_controls_passed':3,'maximum_broad_coverage_identities':{k:v['maximum_reported_identity_among_hits_with90plus_query_coverage'] for k,v in out['databases'].items()},'database_record_count':sum(r['records'] for r in main['indexes']),'all_queries_complete':True,'scope':out['scope']}
 (S/'reference/blast_coverage_summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
