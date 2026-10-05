"""Reconcile the updated panel against completed searches; search missing queries.
The main 133-query result files remain immutable. Merged summaries reference
both successful runs and retain each exact-reference positive control.
"""
from pathlib import Path
import json,hashlib,subprocess,time,sys
from novelty_superset import readseq
S=Path('/mnt/data/egfr_campaign/stage6');B=S/'intermediate/blast';BIN=S/'runtime/bin'
def read(p):return json.loads(Path(p).read_text())
def main():
 t=time.time();original=read(S/'reference/novelty_query_superset.json');seqs={r['sequence'] for r in original};panel=read(S/'reference/panel_preview.json')[:105]
 new=[{'query_id':'Q'+hashlib.sha256(r['sequence'].encode()).hexdigest()[:20],'sequence':r['sequence'],'sources':[r['design_source']]} for r in panel if r['sequence'] not in seqs]
 (S/'reference/novelty_supplement_queries.json').write_text(json.dumps(new,indent=2));assert len({r['sequence'] for r in new})==len(new)
 mergedrows=original+new;(S/'reference/novelty_query_superset_final.json').write_text(json.dumps(mergedrows,indent=2));print('Missing-query supplement',len(new),flush=True)
 for name in ['swissprot','pdb','antibodies_augmented']:
  old=read(S/'reference'/f'panel_blast_{name}_summary.json');summary={'source_main':str(S/'reference'/f'panel_blast_{name}_summary.json'),'database':name,'candidate_queries':len(mergedrows),'all_query_ids_verified':old['all_query_ids_verified'],'positive_control_exact_match':old['positive_control_exact_match'],'hits':dict(old['hits'])}
  if new:
   control=dict(readseq(B/f'panel_superset_{name}.fasta'))['POSITIVE_CONTROL'];records={r['query_id']:r['sequence'] for r in new};records['POSITIVE_CONTROL']=control
   q=B/f'panel_supplement_{name}.fasta';q.write_text(''.join('>'+n+'\n'+seq+'\n' for n,seq in records.items()));dest=B/f'panel_supplement_{name}.tsv'
   cmd=[str(BIN/'blastp'),'-task','blastp','-query',str(q),'-db',str(B/name),'-out',str(dest),'-outfmt','7 qseqid sseqid pident length qlen slen qcovhsp evalue bitscore qstart qend sstart send','-max_target_seqs','50','-num_threads','2','-seg','no','-comp_based_stats','0','-word_size','3','-matrix','BLOSUM62','-evalue','0.001']
   st=time.time();p=subprocess.run(cmd,text=True,capture_output=True,timeout=600);(S/'logs'/f'panel_supplement_{name}_stderr.log').write_text(p.stderr);assert p.returncode==0,p.stderr
   observed=set();hits={qid:[] for qid in records}
   for line in dest.read_text().splitlines():
    if line.startswith('# Query: '):observed.add(line.split(': ',1)[1].split()[0])
    elif line and not line.startswith('#'):
     a=line.split('\t');hits[a[0]].append({'subject':a[1],'identity_pct':float(a[2]),'alignment_length':int(a[3]),'query_length':int(a[4]),'subject_length':int(a[5]),'query_coverage_pct':float(a[6]),'evalue':float(a[7]),'bitscore':float(a[8])})
   assert observed==set(records);cp=any(h['identity_pct']==100 and h['query_coverage_pct']==100 for h in hits['POSITIVE_CONTROL']);assert cp
   extra={'database':name,'candidate_queries':len(new),'all_query_ids_verified':True,'positive_control_exact_match':cp,'command':cmd,'seconds':time.time()-st,'query_sha256':hashlib.sha256(q.read_bytes()).hexdigest(),'result_sha256':hashlib.sha256(dest.read_bytes()).hexdigest(),'hits':hits}
   xp=S/'reference'/f'panel_blast_{name}_supplement_summary.json';xp.write_text(json.dumps(extra,indent=2));summary['source_supplement']=str(xp);summary['hits'].update({k:v for k,v in hits.items() if k!='POSITIVE_CONTROL'})
  assert all(r['query_id'] in summary['hits'] for r in mergedrows)
  (S/'reference'/f'panel_blast_{name}_final_summary.json').write_text(json.dumps(summary,indent=2));print('RECONCILED',name,len(mergedrows),flush=True)
 (S/'reference/panel_blast_final_complete.json').write_text(json.dumps({'complete':True,'main_queries':len(original),'supplement_queries':len(new),'total_candidate_queries':len(mergedrows),'seconds':time.time()-t},indent=2));print('SUPPLEMENT FINISHED',flush=True)
if __name__=='__main__':main()
