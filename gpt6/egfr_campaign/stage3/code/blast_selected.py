"""Whole-sequence BLAST audit of physically evaluated and focused candidates.
Every candidate has separate all-reference CDR3 screening; this is additional
full-protein coverage, not an alternative to the antibody-loop screen.
"""
from pathlib import Path
import json,subprocess,time,hashlib
R=Path('/mnt/data/egfr_campaign');S=R/'stage3'

def main():
    ids={p.name.split('_human6ARU')[0] for p in (S/'intermediate/refined').glob('*_human6ARU_audit.json')}
    ids.update(p.stem for p in (S/'intermediate/designs').glob('C*.json'))
    ids=sorted(ids);out=S/'intermediate/novelty';fasta=out/'physically_selected.fasta';ds=[json.loads((S/'intermediate/designs'/(cid+'.json')).read_text()) for cid in ids]
    # Count sequence-unique entries independently of candidate/model IDs.
    fasta.write_text(''.join('>'+d['candidate_id']+'\n'+d['sequence']+'\n' for d in ds));record={'selection_reason':'All IDs with human MM audit at snapshot, plus every charge-balance design. No selection using novelty results.','ids':ids,'records':len(ds),'unique_sequences':len({d['sequence'] for d in ds}),'query_sha256':hashlib.sha256(fasta.read_bytes()).hexdigest(),'databases':[]};(S/'reference/blast_selected_status.json').write_text(json.dumps(record,indent=2))
    for db in ['swissprot','pdb','antibodies_augmented']:
        target=out/(db+'_selected_hits.tsv');cmd=[str(R/'stage2/runtime/bin/blastp'),'-query',str(fasta),'-db',str(R/'stage2/intermediate/blast'/db),'-out',str(target),'-outfmt','7 qseqid sseqid pident length qlen slen qcovhsp evalue bitscore qstart qend sstart send','-max_target_seqs','10','-num_threads','2','-seg','no','-evalue','0.001'];tic=time.time()
        with open(S/'logs'/('blast_selected_'+db+'.log'),'w') as log:p=subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT)
        row={'database':db,'command':cmd,'exit_code':p.returncode,'seconds':time.time()-tic,'result_sha256':hashlib.sha256(target.read_bytes()).hexdigest() if target.exists() else None};record['databases'].append(row);(S/'reference/blast_selected_status.json').write_text(json.dumps(record,indent=2));print(db,p.returncode,round(row['seconds'],1),flush=True)
        if p.returncode:raise RuntimeError('BLAST failed')
if __name__=='__main__':main()
