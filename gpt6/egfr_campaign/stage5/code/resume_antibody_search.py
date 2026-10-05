"""Resume immutable BLAST query blocks after a controlled runtime interruption.
Scientific search settings are unchanged. Completed outfmt7 query blocks are
reused verbatim. The last open block is always rerun. Positive control required.
"""
from pathlib import Path
import os,signal,time,json,subprocess,sys,hashlib,re
S=Path('/mnt/data/egfr_campaign/stage5');O=S/'intermediate/blast'
def fasta(path):
 out={};key=None
 for l in path.read_text().splitlines():
  if l.startswith('>'):key=l[1:].split()[0];out[key]=''
  else:out[key]+=l.strip()
 return out

def prepare():
 source=O/'antibodies_augmented_stage5.tsv';query=O/'stage5_antibodies_augmented_queries.fasta';records=fasta(query)
 interrupted=[]
 for pid in [1715,1714]:
  try:
   cmd=Path(f'/proc/{pid}/cmdline').read_bytes().decode().replace('\0',' ')
   if 'antibodies_augmented' not in cmd:raise ValueError('PID no longer matches planned worker')
   os.kill(pid,signal.SIGTERM);interrupted.append({'pid':pid,'command':cmd})
  except FileNotFoundError:pass
 time.sleep(1)
 partial=O/'antibodies_augmented_stage5_partial_v1.tsv';source.rename(partial);text=partial.read_text();blocks=[]
 # Each query block begins with a BLASTP comment; preserve only closed blocks.
 spans=[m.start() for m in re.finditer(r'^# BLASTP ',text,re.M)]
 for a,b in zip(spans,spans[1:]):
  block=text[a:b];m=re.search(r'^# Query: (\S+)',block,re.M)
  if m:blocks.append((m.group(1),block))
 done={x[0] for x in blocks};remaining=[k for k in records if k not in done];assert done<=set(records)
 kept=O/'antibodies_completed_blocks.tsv';kept.write_text(''.join(b for k,b in blocks));tasks=[]
 for i in range(0,len(remaining),20):
  ids=remaining[i:i+20];q=O/f'antibody_resume_{i//20:02d}.fasta';q.write_text(''.join('>'+k+'\n'+records[k]+'\n' for k in ids));out=O/f'antibody_resume_{i//20:02d}.tsv';cmd=[str(S/'runtime/bin/blastp'),'-task','blastp','-query',str(q),'-db',str(O/'antibodies_augmented'),'-out',str(out),'-outfmt','7 qseqid sseqid pident length qlen slen qcovhsp evalue bitscore qstart qend sstart send','-max_target_seqs','50','-num_threads','1','-seg','no','-comp_based_stats','0','-word_size','3','-matrix','BLOSUM62','-evalue','0.001'];tasks.append({'query_ids':ids,'query':str(q),'output':str(out),'command':cmd,'query_sha256':hashlib.sha256(q.read_bytes()).hexdigest()})
 plan={'reason':'Long highly homologous antibody search split into finite resumable query chunks; algorithm and search settings unchanged','controlled_interruption':interrupted,'original_query_sha256':hashlib.sha256(query.read_bytes()).hexdigest(),'partial_source':str(partial),'partial_sha256':hashlib.sha256(partial.read_bytes()).hexdigest(),'completed_query_ids':sorted(done),'remaining_query_ids':remaining,'reused_closed_blocks_only':True,'tasks':tasks};(S/'reference/antibody_resume_plan.json').write_text(json.dumps(plan,indent=2));return plan

def execute(plan):
 todo=list(plan['tasks']);active=[];results=[]
 while todo or active:
  while todo and len(active)<2:
   t=todo.pop(0);log=S/'logs'/(Path(t['query']).stem+'.log');f=log.open('w');p=subprocess.Popen(t['command'],stdout=f,stderr=subprocess.STDOUT);active.append((p,f,t,time.time()));print('START',Path(t['query']).name,p.pid,flush=True)
  for item in list(active):
   p,f,t,start=item;rc=p.poll()
   if rc is not None:
    f.close();active.remove(item);results.append({**t,'returncode':rc,'seconds':time.time()-start});print('END',Path(t['query']).name,rc,round(time.time()-start,1),flush=True)
  if active:time.sleep(1)
 assert all(r['returncode']==0 for r in results)
 text=(O/'antibodies_completed_blocks.tsv').read_text()+''.join(Path(t['output']).read_text() for t in plan['tasks']);final=O/'antibodies_augmented_stage5.tsv';final.write_text(text);observed=[];hits={}
 for l in text.splitlines():
  if l.startswith('# Query: '):q=l.split(': ',1)[1].split()[0];observed.append(q);hits.setdefault(q,[])
  elif l and not l.startswith('#'):
   a=l.split('\t');hits[a[0]].append({'subject':a[1],'identity_pct':float(a[2]),'alignment_length':int(a[3]),'query_length':int(a[4]),'subject_length':int(a[5]),'query_coverage_pct':float(a[6]),'evalue':float(a[7]),'bitscore':float(a[8])})
 expected=fasta(O/'stage5_antibodies_augmented_queries.fasta');assert len(observed)==len(set(observed))==len(expected) and set(observed)==set(expected)
 control=any(h['identity_pct']==100 and h['query_coverage_pct']==100 for h in hits['POSITIVE_CONTROL']);assert control
 summary={'database':'antibodies_augmented','query_candidates':len(expected)-1,'observed_query_count_including_control':len(observed),'all_query_ids_verified':True,'positive_control_exact_match':True,'resumed_after_controlled_interruption':True,'plan':str(S/'reference/antibody_resume_plan.json'),'completed_block_count_reused':len(plan['completed_query_ids']),'chunk_results':results,'result_sha256':hashlib.sha256(final.read_bytes()).hexdigest(),'hits':hits,'scope':'Top50 BLAST hits under unchanged original search settings; not exhaustive global alignment.'};(S/'reference/blast_antibodies_augmented_summary.json').write_text(json.dumps(summary,indent=2));print('COMPLETE',len(expected)-1,'queries + exact control',flush=True)
if __name__=='__main__':
 if '--prepare-only' in sys.argv:print('PREPARED',len(prepare()['remaining_query_ids']),flush=True)
 else:execute(json.loads((S/'reference/antibody_resume_plan.json').read_text()))
