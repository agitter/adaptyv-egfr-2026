"""Classical BLAST searches with positive controls and complete query-ID audit."""
from pathlib import Path
import json,sys,subprocess,time,hashlib
R=Path('/mnt/data/egfr_campaign');S=R/'stage5';OUT=S/'intermediate/blast';BIN=S/'runtime/bin'
def run(name):
 start=time.time();designs=[json.loads(p.read_text()) for p in sorted((S/'intermediate/designs').glob('[RN]*.json'))];records={d['candidate_id']:d['sequence'] for d in designs};assert len(records)==len(designs)
 src=OUT/(name+'.fasta')
 from repair_pdb_control import readseq
 header,seq=next((h,s) for h,s in readseq(src) if 80<=len(s)<=500 and set(s)<=set('ACDEFGHIKLMNPQRSTVWY') and (name!='pdb' or 'mol:protein' in h))
 records['POSITIVE_CONTROL']=seq
 query=OUT/f'stage5_{name}_queries.fasta';query.write_text(''.join('>'+k+'\n'+v+'\n' for k,v in records.items()));dest=OUT/f'{name}_stage5.tsv'
 cmd=[str(BIN/'blastp'),'-task','blastp','-query',str(query),'-db',str(OUT/name),'-out',str(dest),'-outfmt','7 qseqid sseqid pident length qlen slen qcovhsp evalue bitscore qstart qend sstart send','-max_target_seqs','50','-num_threads','1','-seg','no','-comp_based_stats','0','-word_size','3','-matrix','BLOSUM62','-evalue','0.001']
 result=subprocess.run(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True,timeout=1200);(S/'logs'/f'blast_{name}_stderr.log').write_text(result.stderr);assert result.returncode==0,result.stderr
 observed=set();hits={q:[] for q in records}
 for line in dest.read_text().splitlines():
  if line.startswith('# Query: '):observed.add(line.split(': ',1)[1].split()[0])
  elif line and not line.startswith('#'):
   a=line.split('\t');hits[a[0]].append({'subject':a[1],'identity_pct':float(a[2]),'alignment_length':int(a[3]),'query_length':int(a[4]),'subject_length':int(a[5]),'query_coverage_pct':float(a[6]),'evalue':float(a[7]),'bitscore':float(a[8])})
 assert observed==set(records),(observed^set(records))
 control=any(h['identity_pct']==100 and h['query_coverage_pct']==100 for h in hits['POSITIVE_CONTROL']);assert control,'Exact-reference positive control failed'
 summary={'database':name,'query_candidates':len(designs),'observed_query_count_including_control':len(observed),'all_query_ids_verified':True,'positive_control_exact_match':control,'positive_control_source':header,'command':cmd,'query_sha256':hashlib.sha256(query.read_bytes()).hexdigest(),'result_sha256':hashlib.sha256(dest.read_bytes()).hexdigest(),'seconds':time.time()-start,'hits':hits,'scope':'Top50 BLAST hits, not exhaustive global alignments; frameworks expected to match known antibodies.'};(S/'reference'/f'blast_{name}_summary.json').write_text(json.dumps(summary,indent=2));print(name,'COMPLETE',len(designs),'seconds',summary['seconds'],flush=True)
if __name__=='__main__':run(sys.argv[1])
