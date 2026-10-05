#!/usr/bin/env python3
"""Fetch public EGFR-campaign inputs; install or execute nothing.

Run: python3 fetch_egfr_inputs.py
Upload: egfr_inputs.zip

Uses only the Python standard library. Reuses valid downloads in legacy_fetch/
from the superseded fetch_legacy_inputs.py script. No Python 2 or historical
BLAST executable is downloaded. BLAST+ 2.17.0 is pinned for reproducibility;
modern implementations are allowed, but design methods remain pre-2011.
"""
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import sys
import time
from urllib.request import Request, urlopen
import zipfile

ROOT = Path('legacy_fetch').resolve()
LOG = ROOT / 'download_v2.log'
VERSION = '2.17.0'
BLAST = 'ncbi-blast-' + VERSION + '+-x64-linux.tar.gz'
BLAST_BASES = [
    'https://ftp.ncbi.nlm.nih.gov/blast/executables/blast%2B/' + VERSION + '/',
    'https://ftp.ncbi.nlm.nih.gov/blast/executables/blast%2B/LATEST/',
]
ITEMS = [
    (BLAST, [b + BLAST.replace('+', '%2B') for b in BLAST_BASES],
     'Modern BLAST implementation for conventional sequence-similarity screening; Linux x86-64 container, not the downloader computer', True),
    (BLAST + '.md5', [b + BLAST.replace('+', '%2B') + '.md5' for b in BLAST_BASES],
     'Publisher checksum for the BLAST archive', False),
    ('uniprot_sprot.fasta.gz', [
        'https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz',
        'https://ftp.expasy.org/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz'],
     'Current reviewed Swiss-Prot sequences for novelty screening, not design seeds', True),
    ('pdb_seqres.txt.gz', [
        'https://files.wwpdb.org/pub/pdb/derived_data/pdb_seqres.txt.gz',
        'https://files.rcsb.org/pub/pdb/derived_data/pdb_seqres.txt.gz'],
     'Current PDB sequences for novelty screening and calibration', True),
    ('uniprot_release_readme.txt', [
        'https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/README'],
     'Swiss-Prot dataset documentation', False),
    ('uniprot_reldate.txt', [
        'https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/reldate.txt'],
     'Swiss-Prot dataset release identifier and date', False),
]
for pdb, role in [
    ('1IVO', '2002 EGFR target conformation; ligand is not a design seed'),
    ('1NQL', '2003 EGFR target conformation; ligand is not a design seed'),
    ('1YY9', '2005 EGFR target conformation; antibody is an excluded reference/control'),
    ('2A3D', '1999 designed-fold calibration control, not a sequence seed'),
    ('1FSD', '1997 designed-fold calibration control, not a sequence seed'),
    ('1L2Y', '2002 designed-fold calibration control, not a sequence seed'),
    ('3DWT', '2008 generic VHH framework candidate; original CDRs must be discarded if used'),
    ('3EAK', '2008 humanized VHH framework candidate; original CDRs must be discarded if used'),
    ('3EBA', '2008 humanized VHH framework candidate; original CDRs must be discarded if used'),
]:
    ITEMS.append((pdb + '.pdb', [
        'https://files.rcsb.org/download/' + pdb + '.pdb',
        'https://www.ebi.ac.uk/pdbe/entry-files/download/pdb' + pdb.lower() + '.ent'],
        role, False))


def utc():
    return time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())


def message(text):
    print(text, flush=True)
    with LOG.open('a', encoding='utf-8') as handle:
        handle.write(utc() + ' ' + text + '\n')


def checksums(path):
    # MD5 is included only to compare NCBI's published checksum, not as a
    # cryptographic authenticity guarantee. SHA-256 is used for provenance.
    sha = hashlib.sha256()
    md5 = hashlib.md5()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            sha.update(block)
            md5.update(block)
    return {'sha256': sha.hexdigest(), 'md5': md5.hexdigest()}


def validate(path, name):
    if not path.is_file() or path.stat().st_size < 20:
        raise ValueError('Missing or unexpectedly small file')
    with path.open('rb') as handle:
        head = handle.read(8192)
    if name.endswith('.gz'):
        if not head.startswith(b'\x1f\x8b'):
            raise ValueError('Not a gzip file')
        # Read through the stream to check gzip completeness and CRC, without
        # creating an uncompressed file or extracting executable contents.
        with gzip.open(path, 'rb') as handle:
            first = handle.read(8192)
            if name in ('uniprot_sprot.fasta.gz', 'pdb_seqres.txt.gz'):
                if not first.lstrip().startswith(b'>'):
                    raise ValueError('Compressed data is not FASTA')
            while handle.read(4 * 1024 * 1024):
                pass
    elif name.endswith('.pdb'):
        if not any(line.startswith((b'HEADER', b'ATOM  ')) for line in head.splitlines()):
            raise ValueError('Not a PDB file')
    elif name.endswith('.md5'):
        if not re.match(rb'\s*[0-9a-fA-F]{32}\b', head):
            raise ValueError('Not a publisher MD5 checksum')
    elif b'<html' in head.lower() or b'<!doctype html' in head.lower():
        raise ValueError('Received HTML instead of the requested data')


def fetch(item):
    name, urls, role, required = item
    path = ROOT / name
    sidecar = ROOT / (name + '.provenance.json')
    rec = {'filename': name, 'role': role, 'required': required,
           'requested_urls': urls, 'attempts': []}
    if path.exists():
        try:
            validate(path, name)
            sums = checksums(path)
            previous = None
            if sidecar.exists():
                previous = json.loads(sidecar.read_text(encoding='utf-8'))
                if previous.get('sha256') and previous['sha256'] != sums['sha256']:
                    raise ValueError('Cached file differs from its provenance checksum')
            rec.update(status='cached', bytes=path.stat().st_size,
                       checked_utc=utc(), **sums)
            if previous is not None:
                rec['prior_provenance'] = previous
            else:
                rec['cache_origin'] = 'Unknown: cached file had no provenance sidecar'
            message('Using validated cached file: ' + name)
            sidecar.write_text(json.dumps(rec, indent=2) + '\n', encoding='utf-8')
            return rec
        except Exception as exc:
            rec['attempts'].append({'stage': 'cache_validation', 'error': str(exc)})
            message('Cached file needs replacement: ' + name)
    rec['status'] = 'failed'
    for url in urls:
        for attempt in range(1, 3):
            partial = ROOT / (name + '.part')
            try:
                message('Downloading ' + name + ' (attempt ' + str(attempt) + ')')
                request = Request(url, headers={'User-Agent': 'EGFRInputFetcher/2.0'})
                with urlopen(request, timeout=120) as response:
                    headers = dict(response.headers.items())
                    final_url = response.geturl()
                    total = 0
                    report_at = 50 * 1024 * 1024
                    with partial.open('wb') as handle:
                        while True:
                            block = response.read(1024 * 1024)
                            if not block:
                                break
                            handle.write(block)
                            total += len(block)
                            if total >= report_at:
                                message('  ' + name + ': ' + str(total // (1024 * 1024)) + ' MiB')
                                report_at += 50 * 1024 * 1024
                    length = response.headers.get('Content-Length')
                    if length and length.isdigit() and total != int(length):
                        raise ValueError('Incomplete HTTP response: byte count differs from Content-Length')
                validate(partial, name)
                os.replace(str(partial), str(path))
                rec.update(status='downloaded', source_url=url, final_url=final_url,
                           retrieved_utc=utc(), http_headers=headers,
                           bytes=path.stat().st_size, **checksums(path))
                sidecar.write_text(json.dumps(rec, indent=2) + '\n', encoding='utf-8')
                return rec
            except Exception as exc:
                rec['attempts'].append({'url': url, 'attempt': attempt,
                                        'error': type(exc).__name__ + ': ' + str(exc)})
                message('  Failed: ' + str(exc))
                if partial.exists():
                    partial.unlink()
                time.sleep(1)
    return rec


def main():
    ROOT.mkdir(parents=True, exist_ok=True)
    message('Fetching public data and BLAST+ 2.17.0. Nothing downloaded will be installed or executed.')
    message('The Linux BLAST archive is for the remote computation container, regardless of your computer OS.')
    message('Previously downloaded Python 2 and BLAST 2.2.24 files are ignored and will not be packaged.')
    records = [fetch(item) for item in ITEMS]
    by_name = {rec['filename']: rec for rec in records}
    binary = by_name[BLAST]
    md5rec = by_name[BLAST + '.md5']
    if binary['status'] != 'failed' and md5rec['status'] != 'failed':
        expected = (ROOT / (BLAST + '.md5')).read_text(encoding='utf-8').split()[0].lower()
        binary['publisher_md5_match'] = binary['md5'] == expected
        if not binary['publisher_md5_match']:
            binary['status'] = 'failed'
            binary['error'] = 'Publisher checksum mismatch; archive excluded'
            message('WARNING: BLAST checksum mismatch; excluding the binary archive.')
        else:
            message('BLAST publisher checksum verified.')
    elif binary['status'] != 'failed':
        binary['publisher_md5_match'] = None
        message('WARNING: Publisher MD5 unavailable; SHA-256 and gzip integrity are recorded, but publisher checksum was not verified.')
    for rec in records:
        (ROOT / (rec['filename'] + '.provenance.json')).write_text(
            json.dumps(rec, indent=2) + '\n', encoding='utf-8')
    manifest_path = ROOT / 'manifest_v2.json'
    manifest_path.write_text(json.dumps({
        'created_utc': utc(), 'fetcher_version': 2, 'files': records,
        'scope': 'Modern implementation allowed; scientific methods must predate 2011.',
        'note': 'No downloaded software was installed or executed; no archive was extracted.',
    }, indent=2) + '\n', encoding='utf-8')
    archive_path = Path('egfr_inputs.zip').resolve()
    message('Packaging ' + str(archive_path))
    with zipfile.ZipFile(archive_path, 'w', zipfile.ZIP_DEFLATED, allowZip64=True) as archive:
        for rec in records:
            name = rec['filename']
            if rec['status'] != 'failed':
                archive.write(ROOT / name, 'legacy_fetch/' + name,
                              compress_type=zipfile.ZIP_STORED if name.endswith('.gz') else zipfile.ZIP_DEFLATED)
            archive.write(ROOT / (name + '.provenance.json'), 'legacy_fetch/' + name + '.provenance.json')
        archive.write(manifest_path, 'legacy_fetch/manifest_v2.json')
        archive.write(LOG, 'legacy_fetch/download_v2.log')
        archive.write(Path(__file__).resolve(), 'fetch_egfr_inputs.py')
    missing = [rec['filename'] for rec in records if rec['status'] == 'failed' and rec['required']]
    message('DONE: ' + str(archive_path))
    message('Archive size: %.1f MiB' % (archive_path.stat().st_size / (1024 * 1024)))
    if missing:
        message('Important files unavailable: ' + ', '.join(missing))
        message('Upload egfr_inputs.zip anyway; successful downloads and all failure records are included.')
    else:
        message('Upload egfr_inputs.zip to the conversation.')
    return 0


if __name__ == '__main__':
    sys.exit(main())
