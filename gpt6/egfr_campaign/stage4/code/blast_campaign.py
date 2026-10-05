"""Rebuild offline novelty indexes and screen this stage's full sequences.
Classical gapped BLAST1997, BLOSUM621992. Version4 indexes; modern executable
per the user's explicit implementation-version allowance. No learned models.
"""
from pathlib import Path
import sys,os,gzip,csv,json,re,hashlib,tarfile,shutil,subprocess,time
R=Path('/mnt/data/egfr_campaign');S=R/'stage4';OUT=S/'intermediate/blast';OUT.mkdir(parents=True,exist_ok=True);BIN=S/'runtime/bin';BIN.mkdir(parents=True,exist_ok=True);D=R/'inputs/transfer2/egfr_refinement_fetch';T=R/'inputs/transfer1/legacy_fetch';AA=set('ACDEFGHIKLMNPQRSTVWYXBZ')
csv.field_size_limit(32*1024*1024)
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
 return h.hexdigest()
def run(cmd,name):
 start=time.time();log=S/'logs'/(name+'.log')
 with log.open('w') as f:r=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,timeout=1200)
 row={'command':list(map(str,cmd)),'returncode':r.returncode,'seconds':time.time()-start,'log':str(log)}
 if r.returncode:raise RuntimeError(str(row))
 return row

def main():
 start=time.time();record={'method':'Gapped BLAST1997, BLOSUM621992; composition adjustment disabled; modern BLAST2.17.0 implementation allowed','indexes':[],'searches':[]};status=S/'reference/blast_stage4_status.json'
 archive=T/'ncbi-blast-2.17.0+-x64-linux.tar.gz'
 with tarfile.open(archive,'r:gz') as t:
  for m in t:
   name=Path(m.name).name
   if name not in ['blastp','makeblastdb'] or not m.name.endswith('/bin/'+name):continue
   assert m.isfile();dest=BIN/name
   if not dest.exists():
    with t.extractfile(m) as src,dest.open('wb') as f:shutil.copyfileobj(src,f)
    dest.chmod(0o755)
 record['binary_archive_sha256']=sha(archive);record['executable_sha256']={n:sha(BIN/n) for n in ['blastp','makeblastdb']}
 record['version']=subprocess.check_output([str(BIN/'blastp'),'-version'],text=True);status.write_text(json.dumps(record,indent=2))
 for name,fname in [('swissprot','uniprot_sprot.fasta.gz'),('pdb','pdb_seqres.txt.gz')]:
  out=OUT/(name+'.fasta');count=0
  if not out.exists():
   with gzip.open(T/fname,'rt') as src,out.open('w') as f:
    for line in src:
     # Includes all supplied PDB polymer records, matching prior conservative
     # index scope; nonprotein records can add false positives, not bypasses.
     if line.startswith('>'):count+=1
     f.write(line)
  else:
   with out.open() as f:count=sum(line.startswith('>') for line in f)
  record['indexes'].append({'database':name,'records':count,'source':str(T/fname),'source_sha256':sha(T/fname),'normalized_fasta_sha256':sha(out)});print('Prepared',name,count,flush=True)
 ab=OUT/'antibodies_augmented.fasta';sources=OUT/'antibody_sources.jsonl';seqs={}
 if not ab.exists():
  def add(seq,src):
   seq=re.sub(r'\s','',str(seq)).upper()
   if len(seq)<30 or not set(seq)<=AA:return
   if seq not in seqs:seqs[seq]=src
  counts={}
  for fname in ['plabdab_paired.csv.gz','plabdab_unpaired.csv.gz','plabdab_nano.csv.gz','therasabdab.csv']:
   f=gzip.open(D/fname,'rt') if fname.endswith('.gz') else open(D/fname,encoding='utf-8-sig');n=0
   for n,row in enumerate(csv.DictReader(f),1):
    src=f'{fname}:row{n}'
    if fname=='plabdab_paired.csv.gz':
     for k in ['heavy_sequence','light_sequence']:add(row.get(k,''),src+':'+k)
    elif fname=='plabdab_unpaired.csv.gz':
     add(row.get('GBSeq_sequence',''),src);add(row.get('numbered','').replace('-',''),src+':provided_numbered_domain')
    elif fname=='plabdab_nano.csv.gz':add(row.get('sequence',''),src)
    else:
     for k in ['HeavySequence','LightSequence','HeavySequence(ifbispec)','LightSequence(ifbispec)']:add(row.get(k,''),src+':'+k)
   f.close();counts[fname]=n;print('Antibody inputs',fname,n,len(seqs),flush=True)
  with ab.open('w') as f,sources.open('w') as g:
   for seq,src in sorted(seqs.items()):
    ident=hashlib.sha256(seq.encode()).hexdigest();f.write('>'+ident+'\n'+seq+'\n');g.write(json.dumps({'id':ident,'first_source':src})+'\n')
  record['antibody_input_rows']=counts;count=len(seqs);seqs.clear()
 else:
  with ab.open() as f:count=sum(line.startswith('>') for line in f)
 record['indexes'].append({'database':'antibodies_augmented','records':count,'normalized_fasta_sha256':sha(ab),'source_map_sha256':sha(sources)});status.write_text(json.dumps(record,indent=2))
 query=OUT/'stage4_INTERNAL_NOT_FOR_SUBMISSION.fasta';designs=[json.loads(p.read_text()) for p in sorted((S/'intermediate/designs').glob('*.json'))];query.write_text(''.join('>'+d['candidate_id']+'\n'+d['sequence']+'\n' for d in designs));record['query_count']=len(designs);record['query_sha256']=sha(query)
 for item in record['indexes']:
  name=item['database'];prefix=OUT/name
  if not prefix.with_suffix('.pin').exists():item['index_build']=run([str(BIN/'makeblastdb'),'-in',str(OUT/(name+'.fasta')),'-out',str(prefix),'-dbtype','prot','-blastdb_version','4'],'build_'+name)
  hits=OUT/(name+'_stage4.tsv');cmd=[str(BIN/'blastp'),'-task','blastp','-query',str(query),'-db',str(prefix),'-out',str(hits),'-outfmt','7 qseqid sseqid pident length qlen slen qcovhsp evalue bitscore qstart qend sstart send','-max_target_seqs','20','-num_threads','1','-seg','no','-comp_based_stats','0','-word_size','3','-matrix','BLOSUM62','-evalue','0.001'];row=run(cmd,'blast_stage4_'+name);row.update(database=name,result_sha256=sha(hits));record['searches'].append(row);status.write_text(json.dumps(record,indent=2));print('Completed',name,round(row['seconds'],1),flush=True)
 record['seconds']=time.time()-start;record['complete']=True;record['scope']='Whole-chain similarity search only; generic antibody-framework hits are expected, and CDR novelty is evaluated separately. Not official Proteinbase classifier or patent-wide eligibility.';status.write_text(json.dumps(record,indent=2));print('ALL DONE',record['seconds'],flush=True)
if __name__=='__main__':main()
