"""Classical Levenshtein (1965) screening and BLAST (1990/1997) preparation.
Uses a modern compiled implementation of the historical edit-distance metric,
as expressly permitted by the user. No learned antibody-numbering tool is run.
"""
from pathlib import Path
import json,sys,time,math,re,subprocess
from rapidfuzz.distance import Levenshtein
from rapidfuzz import process
R=Path('/mnt/data/egfr_campaign');S=R/'stage3'
D=S/'intermediate/designs';E=S/'intermediate/evaluation';N=R/'stage2/intermediate/novelty';OUT=S/'intermediate/novelty';OUT.mkdir(exist_ok=True)

def main():
    refs=[json.loads(line) for line in open(N/'cdrh3_augmented.jsonl')]
    seqs=[r['sequence'] for r in refs];lengths=sorted(set(map(len,seqs)))
    bylen={n:[r for r in refs if len(r['sequence'])==n] for n in lengths};seqby={n:[r['sequence'] for r in bylen[n]] for n in lengths}
    ds=[json.loads(p.read_text()) for p in sorted(D.glob('*.json')) if '_summary' not in p.name and '_failed' not in p.name]
    cache={};out=[];t0=time.time()
    for d in ds:
        seq=d['cdr_sequences'][2];n=len(seq)
        if seq not in cache:
            best=-1.;bestref=None;distbest=None
            for m in lengths:
                # Upper bound on possible identity for different sequence lengths.
                if min(m,n)/max(m,n)<=best:continue
                hit=process.extractOne(seq,seqby[m],scorer=Levenshtein.distance)
                if hit is None:continue
                hitseq,dist,idx=hit;identity=1.-dist/max(n,m)
                if identity>best:best=identity;bestref=bylen[m][idx];distbest=dist
            cache[seq]={'cdr3_query':seq,'closest_cdr3':bestref,'edit_distance':int(distbest),'max_edit_identity':best,'local_cdr3_identity_below_70pct':best<.70}
        out.append({'candidate_id':d['candidate_id'],**cache[seq]})
    (OUT/'candidate_cdr3_screen_augmented.json').write_text(json.dumps(out,indent=2))
    with open(OUT/'candidates.fasta','w') as f:
        for d in ds:f.write('>'+d['candidate_id']+'\n'+d['sequence']+'\n')
    # Reference-positive and near-match controls; independently checked DP.
    def dp(a,b):
        old=list(range(len(b)+1))
        for i,c in enumerate(a,1):
            cur=[i]
            for j,e in enumerate(b,1):cur.append(min(cur[-1]+1,old[j]+1,old[j-1]+(c!=e)))
            old=cur
        return old[-1]
    tests=[]
    for a,b in [('ARDYYGMDV','ARDYYGMDV'),('ARDYYGMDV','ARDYFGMDV'),('AAAA','YYYY'),('HGDHGEYY','DHGEY'),('','AAAA')]:
        assert dp(a,b)==Levenshtein.distance(a,b);tests.append([a,b,dp(a,b)])
    probe=seqs[len(seqs)//2];assert process.extractOne(probe,seqs,scorer=Levenshtein.distance)[1]==0
    report={'candidates':len(ds),'unique_cdr3_queries':len(cache),'reference_cdr3':len(refs),'below_70pct':sum(a['local_cdr3_identity_below_70pct'] for a in out),'max_identity':max(a['max_edit_identity'] for a in out),'distance_tests':tests,'positive_control_exact_match':True,'seconds':time.time()-t0,'qualification':'Local, motif-based coverage only; not organizer eligibility certification.'}
    (S/'reference/novelty_screen_augmented_stats.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
if __name__=='__main__':main()
