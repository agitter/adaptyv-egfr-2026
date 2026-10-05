"""Additional query coverage plus exact database positive controls."""
from pathlib import Path
import json,hashlib,subprocess,time
R=Path('/mnt/data/egfr_campaign');S=R/'stage4';OUT=S/'intermediate/blast';AA=set('ACDEFGHIKLMNPQRSTVWY')
def first_control(p):
 ident=None;seq=''
 with open(p) as f:
  for line in f:
   if line.startswith('>'):
    if ident is not None and 40<=len(seq)<=200 and set(seq)<=AA:return ident,seq
    ident=line[1:].strip().split()[0];seq=''
   else:seq+=line.strip()
 raise ValueError('No appropriate positive-control sequence')
def main():
 out=[]
 for name in ['swissprot','pdb','antibodies_augmented']:
  ident,seq=first_control(OUT/(name+'.fasta'));query=OUT/(name+'_followup_queries.fasta');ids=['S00062','S00063'];text=''
  for cid in ids:
   d=json.loads((S/'intermediate/designs'/(cid+'.json')).read_text());text+='>'+cid+'\n'+d['sequence']+'\n'
  control='SELF_CHECK_'+name;text+='>'+control+'\n'+seq+'\n';query.write_text(text);hits=OUT/(name+'_followup.tsv');cmd=[str(S/'runtime/bin/blastp'),'-task','blastp','-query',str(query),'-db',str(OUT/name),'-out',str(hits),'-outfmt','7 qseqid sseqid pident length qlen slen qcovhsp evalue bitscore qstart qend sstart send','-max_target_seqs','20','-num_threads','1','-seg','no','-comp_based_stats','0','-word_size','3','-matrix','BLOSUM62','-evalue','0.001'];start=time.time()
  with open(S/'logs'/('blast_followup_'+name+'.log'),'w') as f:p=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,timeout=600)
  assert p.returncode==0;rows=[line.split('\t') for line in hits.read_text().splitlines() if line and not line.startswith('#')];positive=any(r[0]==control and float(r[2])==100. and int(r[3])==len(seq) and int(r[4])==len(seq) for r in rows);assert positive,name
  record={'database':name,'command':cmd,'returncode':p.returncode,'seconds':time.time()-start,'additional_candidate_ids':ids,'positive_control_query_id':control,'positive_control_reference_id':ident,'positive_control_length':len(seq),'positive_control_exact_full_match_pass':positive,'query_sha256':hashlib.sha256(query.read_bytes()).hexdigest(),'results_sha256':hashlib.sha256(hits.read_bytes()).hexdigest()};out.append(record);(S/'reference/blast_followup_status.json').write_text(json.dumps(out,indent=2));print(name,'positive control',positive,'seconds',record['seconds'],flush=True)
if __name__=='__main__':main()
