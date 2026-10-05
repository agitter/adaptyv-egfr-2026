"""Retain failed DNA-record control and run a valid protein control separately.
All 205 original sequence searches are preserved; no redundant full rerun.
"""
from pathlib import Path
import json,subprocess,hashlib,time
S=Path('/mnt/data/egfr_campaign/stage5');O=S/'intermediate/blast'
def readseq(path):
 header=None;seq=[]
 with path.open() as f:
  for l in f:
   if l.startswith('>'):
    if header is not None:yield header,''.join(seq)
    header=l[1:].strip();seq=[]
   else:seq.append(l.strip())
  if header is not None:yield header,''.join(seq)
def main():
 header,seq=next((h,s) for h,s in readseq(O/'pdb.fasta') if 'mol:protein' in h and 80<=len(s)<=500 and set(s)<=set('ACDEFGHIKLMNPQRSTVWY'))
 query=O/'pdb_valid_protein_control.fasta';query.write_text('>PDB_PROTEIN_CONTROL\n'+seq+'\n');dest=O/'pdb_valid_protein_control.tsv'
 cmd=[str(S/'runtime/bin/blastp'),'-task','blastp','-query',str(query),'-db',str(O/'pdb'),'-out',str(dest),'-outfmt','7 qseqid sseqid pident length qlen slen qcovhsp evalue bitscore qstart qend sstart send','-max_target_seqs','50','-num_threads','1','-seg','no','-comp_based_stats','0','-word_size','3','-matrix','BLOSUM62','-evalue','0.001']
 t=time.time();p=subprocess.run(cmd,text=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,timeout=180);(S/'logs/pdb_corrected_control_stderr.log').write_text(p.stderr);assert p.returncode==0
 control=[]
 for l in dest.read_text().splitlines():
  if l and not l.startswith('#'):control.append(l.split('\t'))
 assert any(float(x[2])==100. and float(x[6])==100. for x in control)
 observed=set();hits={}
 for l in (O/'pdb_stage5.tsv').read_text().splitlines():
  if l.startswith('# Query: '):q=l.split(': ',1)[1].split()[0];observed.add(q);hits.setdefault(q,[])
  elif l and not l.startswith('#'):
   a=l.split('\t');hits[a[0]].append({'subject':a[1],'identity_pct':float(a[2]),'alignment_length':int(a[3]),'query_length':int(a[4]),'subject_length':int(a[5]),'query_coverage_pct':float(a[6]),'evalue':float(a[7]),'bitscore':float(a[8])})
 expected={p.stem for p in (S/'intermediate/designs').glob('[RN]*.json')};assert observed==expected|{'POSITIVE_CONTROL'}
 r={'database':'pdb','query_candidates':len(expected),'observed_query_count_including_control':len(observed),'all_query_ids_verified':True,'positive_control_exact_match':True,'positive_control_source':header,'command':cmd,'control_seconds':time.time()-t,'hits':hits,'original_control_failure':{'source':'100d_A mol:na length:10','reason':'First PDB record was a 10-base nucleic-acid polymer, inappropriate protein-search control','preserved_original_output':str(O/'pdb_stage5.tsv'),'original_control_exact_match':False},'corrected_control_hits':control,'files':[{ 'path':str(q),'sha256':hashlib.sha256(q.read_bytes()).hexdigest()} for q in [query,dest,O/'pdb_stage5.tsv',O/'stage5_pdb_queries.fasta']],'scope':'205 original query searches verified; valid protein positive control added without replacing original failed control.'}
 (S/'reference/blast_pdb_summary.json').write_text(json.dumps(r,indent=2));print('Corrected protein control PASS',header,'205 query IDs verified',flush=True)
if __name__=='__main__':main()
