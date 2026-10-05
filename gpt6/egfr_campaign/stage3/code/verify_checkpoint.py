"""Independently hash every archived payload against the archived manifest."""
from pathlib import Path
import zipfile,json,hashlib,datetime
p=Path('/mnt/data/egfr_stage3_checkpoint.zip');S=Path('/mnt/data/egfr_campaign/stage3')
with zipfile.ZipFile(p) as z:
    rows=json.loads(z.read('egfr_campaign/stage3/reference/checkpoint_included_manifest.json'))
    for row in rows:
        name='egfr_campaign/'+row['path'];h=hashlib.sha256();size=0
        with z.open(name) as f:
            for b in iter(lambda:f.read(1024*1024),b''):h.update(b);size+=len(b)
        assert h.hexdigest()==row['sha256'],name
        assert size==row['bytes'],name
    assert len(z.namelist())==len(rows)+3
result={'archive':str(p),'verified_payload_files':len(rows),'all_payload_sha256_and_sizes_match':True,'archive_members':len(rows)+3,'timestamp_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
(S/'reference/checkpoint_independent_verification.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
