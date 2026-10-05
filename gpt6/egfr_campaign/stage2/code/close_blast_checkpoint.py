"""Stop only this checkpoint's obsolete broad-antibody search; record scope."""
from pathlib import Path
import json,os,signal,time,subprocess
S=Path('/mnt/data/egfr_campaign/stage2')
# The completed exhaustive CDR screen already uses the augmented references.
# This broad search uses the INITIAL index and will be replaced after redesign.
ps=subprocess.check_output(['ps','-eo','pid,args'],text=True)
stopped=[]
for line in ps.splitlines()[1:]:
    fields=line.strip().split(None,1)
    if len(fields)!=2:continue
    pid=int(fields[0]);cmd=fields[1]
    relevant=(cmd.startswith(str(S/'runtime/bin/blastp')+' ') and str(S/'intermediate/blast/antibodies') in cmd) or cmd==('sh '+str(S/'code/search_novelty.sh'))
    if relevant:
        try:os.kill(pid,signal.SIGTERM);stopped.append({'pid':pid,'command':cmd})
        except ProcessLookupError:pass
queries={l[1:].split()[0] for l in (S/'intermediate/novelty/candidates.fasta').read_text().splitlines() if l.startswith('>')}
record={'query_count':len(queries),'stopped_processes':stopped,'searches':{},'reason_for_stopping':'Checkpoint boundary; broad antibody job uses superseded initial reference index. Exhaustive augmented CDR3 screen is complete. Repeat broad antibody search against augmented index after sequence redesign.'}
for name in ['swissprot','pdb','antibodies']:
    p=S/'intermediate/novelty'/('blast_'+name+'.tsv')
    lines=p.read_text().splitlines() if p.exists() else []
    q={l.split('\t')[0] for l in lines if l.strip()}
    item={'rows':len(lines),'queries_with_hits':len(q),'reference_count':{'swissprot':575748,'pdb':1165667,'antibodies':374451}[name]}
    if name in ['swissprot','pdb']:
        assert q==queries,(name,len(q))
        item['status']='complete';item['completion_evidence']='sequential set -e driver advanced to later search; all query identifiers present'
        # Full-query coverage matters more than short local high-identity matches.
        hits=[l.split('\t') for l in lines]
        covered=[h for h in hits if float(h[12])>=80.]
        item['maximum_percent_identity_among_hits_with_80pct_query_coverage']=max((float(h[2]) for h in covered),default=None)
        item['interpretation']='Top-20 local BLAST hits; common antibody frameworks are expected. Not a full organizer novelty decision.'
    else:
        item['status']='stopped_incomplete'
        if p.exists():
            partial=p.with_name(p.stem+'.partial.tsv');p.rename(partial);item['partial_output']=str(partial.relative_to(S))
    record['searches'][name]=item
(S/'reference/blast_status.json').write_text(json.dumps(record,indent=2));print(json.dumps(record,indent=2))
