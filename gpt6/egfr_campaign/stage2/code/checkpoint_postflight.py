"""Independently verify completed scientific output consistency and archive hashes."""
from pathlib import Path
import json,math,hashlib,zipfile,ast,sys
S=Path('/mnt/data/egfr_campaign/stage2')

def finite(x):
    if isinstance(x,float):return math.isfinite(x)
    if isinstance(x,dict):return all(finite(v) for v in x.values())
    if isinstance(x,list):return all(finite(v) for v in x)
    return True

def outputs():
    scripts=list((S/'code').glob('*.py'))
    for p in scripts:ast.parse(p.read_text(),filename=str(p))
    seq=[];n=0
    for p in (S/'intermediate/designs').glob('P*.json'):
        if p.stem.endswith(('_summary','_failed')):continue
        d=json.loads(p.read_text());assert finite(d),p
        s=d['sequence'];assert s==''.join(r['aa'] for r in d['structure']),p
        assert set(s)<=set('ACDEFGHIKLMNPQRSTVWY') and 10<=len(s)<=250
        assert d['cdr_sequences']==[s[a:b] for a,b in d['cdr_intervals_zero_based']],p
        seq.append(s);n+=1
    assert n==366 and len(set(seq))==366
    for p in (S/'intermediate/evaluation').glob('*.json'):assert finite(json.loads(p.read_text())),p
    r={'syntax_checked_scripts':len(scripts),'sequence_structure_cdr_consistency_checked':n,'unique_sequences':len(set(seq)),'finite_numeric_outputs':True,'pass':True}
    (S/'reference/checkpoint_postflight.json').write_text(json.dumps(r,indent=2));print(json.dumps(r))

def archive(path):
    with zipfile.ZipFile(path) as z:
        rows=json.loads(z.read('egfr_campaign/stage2/reference/checkpoint_included_manifest.json'))
        for r in rows:
            data=z.read('egfr_campaign/'+r['path']);assert len(data)==r['bytes'];assert hashlib.sha256(data).hexdigest()==r['sha256'],r['path']
        assert z.testzip() is None
    print(json.dumps({'archive_sha256_manifest_checks':len(rows),'pass':True}))
if __name__=='__main__':
    if len(sys.argv)>1:archive(sys.argv[1])
    else:outputs()
