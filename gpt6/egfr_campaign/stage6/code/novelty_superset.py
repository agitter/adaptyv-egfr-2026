"""Search a final-panel superset with classical BLAST and exact-reference controls.
Antibody whole-chain matches are reported, not treated as a blanket novelty fail.
No organizer classifier, learned predictor, MMseqs2 or modern numbering is run.
"""
from pathlib import Path
import json,hashlib,subprocess,time,sys,os
S=Path('/mnt/data/egfr_campaign/stage6');OUT=S/'intermediate/blast';BIN=S/'runtime/bin';AA=set('ACDEFGHIKLMNPQRSTVWY')
def readseq(p):
 h=None;s=[]
 with Path(p).open() as f:
  for line in f:
   if line.startswith('>'):
    if h is not None:yield h,''.join(s)
    h=line[1:].strip();s=[]
   else:s.append(line.strip())
 if h is not None:yield h,''.join(s)
def prepare():
 records={}
 for r in json.loads((S/'reference/panel_preview.json').read_text()):records[r['sequence']]={'sources':[r['design_source']]}
 for p in (S/'intermediate/designs').glob('[JK]*.json'):
  d=json.loads(p.read_text());records.setdefault(d['sequence'],{'sources':[]})['sources'].append(str(p))
 rows=[]
 for seq,meta in sorted(records.items()):rows.append({'query_id':'Q'+hashlib.sha256(seq.encode()).hexdigest()[:20],'sequence':seq,**meta})
 assert len({r['query_id'] for r in rows})==len(rows)
 (S/'reference/novelty_query_superset.json').write_text(json.dumps(rows,indent=2));print('Prepared',len(rows),'unique candidate queries',flush=True)

def run(name):
 t=time.time();rows=json.loads((S/'reference/novelty_query_superset.json').read_text());records={r['query_id']:r['sequence'] for r in rows};header,seq=next((h,s) for h,s in readseq(OUT/(name+'.fasta')) if 80<=len(s)<=500 and set(s)<=AA and (name!='pdb' or 'mol:protein' in h));records['POSITIVE_CONTROL']=seq
 query=OUT/f'panel_superset_{name}.fasta';query.write_text(''.join('>'+k+'\n'+v+'\n' for k,v in records.items()));dest=OUT/f'panel_superset_{name}.tsv'
 cmd=[str(BIN/'blastp'),'-task','blastp','-query',str(query),'-db',str(OUT/name),'-out',str(dest),'-outfmt','7 qseqid sseqid pident length qlen slen qcovhsp evalue bitscore qstart qend sstart send','-max_target_seqs','50','-num_threads','2','-seg','no','-comp_based_stats','0','-word_size','3','-matrix','BLOSUM62','-evalue','0.001']
 result=subprocess.run(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True,timeout=2400);(S/'logs'/f'panel_blast_{name}_stderr.log').write_text(result.stderr);assert result.returncode==0,result.stderr
 observed=set();hits={q:[] for q in records}
 for line in dest.read_text().splitlines():
  if line.startswith('# Query: '):observed.add(line.split(': ',1)[1].split()[0])
  elif line and not line.startswith('#'):
   a=line.split('\t');hits[a[0]].append({'subject':a[1],'identity_pct':float(a[2]),'alignment_length':int(a[3]),'query_length':int(a[4]),'subject_length':int(a[5]),'query_coverage_pct':float(a[6]),'evalue':float(a[7]),'bitscore':float(a[8])})
 assert observed==set(records),(len(observed),len(records),observed^set(records));control=any(h['identity_pct']==100 and h['query_coverage_pct']==100 for h in hits['POSITIVE_CONTROL']);assert control
 summary={'database':name,'candidate_queries':len(rows),'observed_queries_including_control':len(observed),'all_query_ids_verified':True,'positive_control_exact_match':control,'positive_control_source':header,'command':cmd,'query_sha256':hashlib.sha256(query.read_bytes()).hexdigest(),'result_sha256':hashlib.sha256(dest.read_bytes()).hexdigest(),'seconds':time.time()-t,'hits':hits,'scope':'Top50 BLAST hits, not exhaustive global alignment. Official antibody CDRH3 novelty is separate; framework familiarity is not rejected.'}
 (S/'reference'/f'panel_blast_{name}_summary.json').write_text(json.dumps(summary,indent=2));print('COMPLETE',name,len(rows),round(summary['seconds'],1),flush=True)
if __name__=='__main__':
 if len(sys.argv)>1 and sys.argv[1]=='prepare':prepare()
 else:
  for name in ['swissprot','pdb','antibodies_augmented']:run(name)
  (S/'reference/panel_blast_complete.json').write_text(json.dumps({'complete':True,'databases':['swissprot','pdb','antibodies_augmented']}));print('ALL PANEL SEARCHES COMPLETE',flush=True)
