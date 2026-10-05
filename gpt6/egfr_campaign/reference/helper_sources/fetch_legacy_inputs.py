#!/usr/bin/env python
"""Fetch public inputs for an offline, pre-2011-method EGFR design campaign.

Run: python3 fetch_legacy_inputs.py
Output: legacy_inputs.zip, legacy_fetch/manifest.json, legacy_fetch/download.log
No external packages, administrator privileges, installation, or execution of
any downloaded software is required. Existing successful downloads are reused.
The language features/standard-library algorithms used here existed by 2010;
the Python 3 compatibility import is only transport infrastructure.
"""
from __future__ import print_function
import os
import sys
import time
import json
import hashlib
import zipfile
try:
    from urllib.request import Request, urlopen
except ImportError:
    from urllib2 import Request, urlopen

ROOT = os.path.abspath('legacy_fetch')
if not os.path.isdir(ROOT):
    os.makedirs(ROOT)
LOG = os.path.join(ROOT, 'download.log')

ITEMS = [
    ('Python-2.7.tgz',
     ['https://www.python.org/ftp/python/2.7/Python-2.7.tgz'],
     'CPython 2.7.0 source, released 2010; optional isolated legacy interpreter', True),
    ('ncbi-blast-2.2.24+-x64-linux.tar.gz',
     ['https://ftp.ncbi.nlm.nih.gov/blast/executables/blast%2B/2.2.24/ncbi-blast-2.2.24%2B-x64-linux.tar.gz',
      'https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.2.24/ncbi-blast-2.2.24+-x64-linux.tar.gz'],
     'NCBI BLAST+ 2.2.24 Linux executables, released 2010; sequence novelty screening', True),
    ('ncbi-blast-2.2.24+-x64-linux.tar.gz.md5',
     ['https://ftp.ncbi.nlm.nih.gov/blast/executables/blast%2B/2.2.24/ncbi-blast-2.2.24%2B-x64-linux.tar.gz.md5'],
     'Publisher checksum for archived BLAST binaries', False),
    ('uniprot_sprot.fasta.gz',
     ['https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz',
      'https://ftp.expasy.org/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz'],
     'Current reviewed Swiss-Prot sequences; novelty checking, never design seeds', True),
    ('pdb_seqres.txt.gz',
     ['https://files.wwpdb.org/pub/pdb/derived_data/pdb_seqres.txt.gz',
      'https://files.rcsb.org/pub/pdb/derived_data/pdb_seqres.txt.gz'],
     'Current PDB sequences; novelty checking and calibration only', True),
    ('uniprot_release_readme.txt',
     ['https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/README'],
     'Dataset release metadata', False),
]
for pdb, role in [
    ('1IVO', '2002 EGFR ligand-bound target conformation; no ligand sequence used as a seed'),
    ('1NQL', '2003 EGFR target conformation; no ligand sequence used as a seed'),
    ('1YY9', '2005 EGFR target conformation; antibody only an excluded reference/control'),
    ('2A3D', '1999 designed-fold calibration control, not a sequence seed'),
    ('1FSD', '1997 designed-fold calibration control, not a sequence seed'),
    ('1L2Y', '2002 designed-fold calibration control, not a sequence seed'),
    ('3DWT', '2008 generic VHH framework candidate; original CDRs must be discarded if used'),
    ('3EAK', '2008 humanized VHH framework candidate; original CDRs must be discarded if used'),
    ('3EBA', '2008 humanized VHH framework candidate; original CDRs must be discarded if used')
]:
    ITEMS.append((pdb + '.pdb',
        ['https://files.rcsb.org/download/' + pdb + '.pdb',
         'https://www.ebi.ac.uk/pdbe/entry-files/download/pdb' + pdb.lower() + '.ent'],
        role, False))


def message(text):
    print(text)
    sys.stdout.flush()
    with open(LOG, 'a') as handle:
        handle.write(time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()) + ' ' + text + '\n')


def digest(path):
    sha = hashlib.sha256()
    with open(path, 'rb') as handle:
        while True:
            block = handle.read(1024 * 1024)
            if not block:
                break
            sha.update(block)
    return sha.hexdigest()


def valid(path, name):
    if not os.path.isfile(path) or os.path.getsize(path) < 20:
        return False
    with open(path, 'rb') as handle:
        head = handle.read(2048)
    if name.endswith(('.gz', '.tgz')):
        return head[:2] == b'\x1f\x8b'
    if name.endswith('.pdb'):
        return b'HEADER' in head or b'ATOM' in head
    return b'<html' not in head.lower() and b'<!doctype html' not in head.lower()


message('Fetching public datasets and archived software. Nothing downloaded will be executed.')
message('Expect a few hundred megabytes of downloads; keep approximately 1 GB of free disk space.')
records = []
for name, urls, role, required in ITEMS:
    path = os.path.join(ROOT, name)
    metadata_path = path + '.provenance.json'
    rec = {'filename': name, 'role': role, 'required': required, 'attempts': []}
    if valid(path, name):
        if os.path.isfile(metadata_path):
            try:
                with open(metadata_path) as handle:
                    rec = json.load(handle)
            except Exception:
                pass
        rec['status'] = 'cached'
        message('Using cached ' + name)
    else:
        rec['status'] = 'failed'
        for url in urls:
            for attempt in range(2):
                partial = path + '.part'
                try:
                    message('Downloading ' + name + ' (attempt ' + str(attempt + 1) + ')')
                    request = Request(url, headers={'User-Agent': 'LegacyInputFetcher/1.0 (public research data)'})
                    response = urlopen(request, timeout=180)
                    headers = dict(response.info().items())
                    downloaded = 0
                    next_report = 25 * 1024 * 1024
                    with open(partial, 'wb') as handle:
                        while True:
                            block = response.read(1024 * 1024)
                            if not block:
                                break
                            handle.write(block)
                            downloaded += len(block)
                            if downloaded >= next_report:
                                message('  ' + name + ': ' + str(downloaded // (1024 * 1024)) + ' MiB')
                                next_report += 25 * 1024 * 1024
                    response.close()
                    if not valid(partial, name):
                        raise ValueError('Response is empty or not the expected file type')
                    if os.path.exists(path):
                        os.remove(path)
                    os.rename(partial, path)
                    rec.update({'status': 'downloaded', 'source_url': url,
                                'retrieved_utc': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
                                'http_headers': headers})
                    break
                except Exception as exc:
                    error = str(exc)
                    rec['attempts'].append({'url': url, 'error': error})
                    message('  Could not fetch: ' + error)
                    if os.path.exists(partial):
                        os.remove(partial)
                    time.sleep(2)
            if rec['status'] == 'downloaded':
                break
    if valid(path, name):
        rec['bytes'] = os.path.getsize(path)
        rec['sha256'] = digest(path)
        with open(metadata_path, 'w') as handle:
            json.dump(rec, handle, indent=2, sort_keys=True)
    records.append(rec)

manifest = {'created_utc': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
            'files': records,
            'note': 'No downloaded code was installed or executed by this script.'}
with open(os.path.join(ROOT, 'manifest.json'), 'w') as handle:
    json.dump(manifest, handle, indent=2, sort_keys=True)
archive_path = os.path.abspath('legacy_inputs.zip')
message('Packaging ' + archive_path)
with zipfile.ZipFile(archive_path, 'w', zipfile.ZIP_DEFLATED, allowZip64=True) as archive:
    for name in sorted(os.listdir(ROOT)):
        path = os.path.join(ROOT, name)
        if os.path.isfile(path) and not name.endswith('.part'):
            method = zipfile.ZIP_STORED if name.endswith(('.gz', '.tgz')) else zipfile.ZIP_DEFLATED
            archive.write(path, 'legacy_fetch/' + name, compress_type=method)
missing = [r['filename'] for r in records if r['status'] == 'failed' and r['required']]
message('DONE: ' + archive_path)
message('Archive size: %.1f MiB' % (os.path.getsize(archive_path) / (1024.0 * 1024.0)))
if missing:
    message('Some important files failed: ' + ', '.join(missing))
    message('Upload the archive anyway; it contains successful files and the failure report.')
else:
    message('Upload legacy_inputs.zip to the conversation.')
