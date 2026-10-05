#!/usr/bin/env python3
"""Phone-a-Friend 2: fetch sequence-only novelty data and classical MM software.

Run with Python 3 (standard library only):
    python fetch_egfr_refinement_inputs.py
Upload the resulting egfr_refinement_inputs.zip.

This downloads but NEVER installs or executes any downloaded software.  The
OpenMM wheel is for the remote Linux x86-64 / CPython 3.13 computation container,
not necessarily for the computer executing this downloader.  Precomputed
antibody models and modern protein-design packages are NOT downloaded.
"""
from __future__ import annotations

import csv
import gzip
import hashlib
from html.parser import HTMLParser
import io
import json
import os
from pathlib import Path
import re
import sys
import time
from urllib.parse import urljoin, urlparse
from urllib.request import Request, urlopen
import zipfile

VERSION = '1.0'
OPENMM_VERSION = '8.4.0.post2'
ROOT = Path('egfr_refinement_fetch').resolve()
LOG = ROOT / 'download.log'
ARCHIVE = Path('egfr_refinement_inputs.zip').resolve()
MAX_DOWNLOAD = 512 * 1024 * 1024
OPENMM_JSON = 'https://pypi.org/pypi/OpenMM/' + OPENMM_VERSION + '/json'
PLABDAB = 'https://opig.stats.ox.ac.uk/webapps/plabdab'
NANO = 'https://opig.stats.ox.ac.uk/webapps/plabdab-nano'


def utc():
    return time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())


def message(text):
    print(text, flush=True)
    ROOT.mkdir(parents=True, exist_ok=True)
    with LOG.open('a', encoding='utf-8') as out:
        out.write(utc() + ' ' + text + '\n')


def sha256(path):
    digest = hashlib.sha256()
    with path.open('rb') as src:
        for block in iter(lambda: src.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def validate(path, filename, expected_sha=None, expected_size=None):
    if not path.is_file() or path.stat().st_size < 20:
        raise ValueError('Missing or unexpectedly small file')
    if expected_size is not None and path.stat().st_size != expected_size:
        raise ValueError('Byte count differs from the publisher metadata')
    digest = sha256(path)
    if expected_sha and digest.lower() != expected_sha.lower():
        raise ValueError('SHA-256 differs from the publisher metadata')
    with path.open('rb') as src:
        head = src.read(8192)
    if filename.endswith('.csv.gz'):
        if not head.startswith(b'\x1f\x8b'):
            raise ValueError('Expected compressed CSV, not gzip data')
        # Validate the entire gzip stream/CRC without extracting a large file.
        with gzip.open(path, 'rb') as src:
            first = src.read(16384)
            if b'<html' in first.lower() or b'<!doctype' in first.lower():
                raise ValueError('Received compressed HTML instead of sequence data')
            first_line = first.decode('utf-8-sig', errors='replace').splitlines()[0]
            header = next(csv.reader([first_line]))
            if len(header) < 2 or not any('seq' in x.lower() or 'heavy' in x.lower() for x in header):
                raise ValueError('CSV header does not look like an antibody sequence dataset')
            total = len(first)
            while True:
                block = src.read(4 * 1024 * 1024)
                if not block:
                    break
                total += len(block)
                if total > 4 * 1024**3:
                    raise ValueError('Unexpectedly large expanded dataset')
    elif filename.endswith('.whl'):
        with zipfile.ZipFile(path) as wheel:
            bad = wheel.testzip()
            if bad:
                raise ValueError('Wheel ZIP checksum failed: ' + bad)
            names = wheel.namelist()
            if not any(x.startswith('openmm/') for x in names):
                raise ValueError('Wheel does not contain the OpenMM package')
            if not any(x.endswith('.dist-info/METADATA') for x in names):
                raise ValueError('Wheel distribution metadata is missing')
    elif filename.endswith('.json'):
        with path.open(encoding='utf-8') as src:
            json.load(src)
    elif filename.endswith('.html'):
        if b'<html' not in head.lower() and b'<!doctype' not in head.lower():
            raise ValueError('Expected an HTML dataset landing page')
    elif filename.endswith('.pdb'):
        if not any(line.startswith((b'HEADER', b'ATOM  ')) for line in head.splitlines()):
            raise ValueError('Expected PDB coordinate data')
    elif b'<html' in head.lower() or b'<!doctype' in head.lower():
        raise ValueError('Received an HTML error page instead of requested data')
    return digest


def fetch(filename, urls, role, required=True, expected_sha=None, expected_size=None):
    if Path(filename).name != filename or filename in ('.', '..'):
        raise ValueError('Unsafe output filename')
    path = ROOT / filename
    sidecar = ROOT / (filename + '.provenance.json')
    rec = dict(filename=filename, requested_urls=list(urls), role=role,
               required=required, attempts=[], status='failed')
    if expected_sha:
        rec['publisher_sha256'] = expected_sha
    if path.exists():
        try:
            digest = validate(path, filename, expected_sha, expected_size)
            prior = json.loads(sidecar.read_text(encoding='utf-8')) if sidecar.exists() else None
            if prior and prior.get('sha256') and prior['sha256'] != digest:
                raise ValueError('Cached file differs from its recorded checksum')
            rec.update(status='cached', bytes=path.stat().st_size,
                       sha256=digest, checked_utc=utc(), prior_provenance=prior)
            message('Using validated cached file: ' + filename)
            sidecar.write_text(json.dumps(rec, indent=2) + '\n', encoding='utf-8')
            return rec
        except Exception as exc:
            rec['attempts'].append(dict(stage='cache_validation', error=str(exc)))
            message('Cached file needs replacement: ' + filename)
    for url in urls:
        for attempt in range(1, 3):
            partial = ROOT / (filename + '.part')
            try:
                if urlparse(url).scheme != 'https':
                    raise ValueError('Only HTTPS source URLs are accepted')
                message('Downloading ' + filename + ' (attempt ' + str(attempt) + ')')
                request = Request(url, headers={'User-Agent': 'EGFRRefinementInputFetcher/' + VERSION})
                with urlopen(request, timeout=120) as response:
                    final_url = response.geturl()
                    if urlparse(final_url).scheme != 'https':
                        raise ValueError('Refusing a redirect to non-HTTPS transport')
                    headers = dict(response.headers.items())
                    total = 0
                    report_at = 25 * 1024 * 1024
                    with partial.open('wb') as out:
                        while True:
                            block = response.read(1024 * 1024)
                            if not block:
                                break
                            total += len(block)
                            if total > MAX_DOWNLOAD:
                                raise ValueError('Download exceeded the per-file size guard')
                            out.write(block)
                            if total >= report_at:
                                message('  ' + filename + ': ' + str(total // (1024 * 1024)) + ' MiB')
                                report_at += 25 * 1024 * 1024
                    length = response.headers.get('Content-Length')
                    if length and length.isdigit() and int(length) != total:
                        raise ValueError('Incomplete HTTP response (Content-Length mismatch)')
                digest = validate(partial, filename, expected_sha, expected_size)
                os.replace(str(partial), str(path))
                rec.update(status='downloaded', source_url=url, final_url=final_url,
                           retrieved_utc=utc(), http_headers=headers,
                           bytes=path.stat().st_size, sha256=digest)
                sidecar.write_text(json.dumps(rec, indent=2) + '\n', encoding='utf-8')
                return rec
            except Exception as exc:
                rec['attempts'].append(dict(url=url, attempt=attempt,
                    error=type(exc).__name__ + ': ' + str(exc)))
                message('  Failed: ' + str(exc))
                if partial.exists():
                    partial.unlink()
                time.sleep(1)
    sidecar.write_text(json.dumps(rec, indent=2) + '\n', encoding='utf-8')
    return rec


class Links(HTMLParser):
    def __init__(self):
        super().__init__(convert_charrefs=True)
        self.links = []
        self.current = None

    def handle_starttag(self, tag, attrs):
        if tag.lower() == 'a':
            self.current = [dict(attrs).get('href', ''), []]

    def handle_data(self, data):
        if self.current is not None:
            self.current[1].append(data)

    def handle_endtag(self, tag):
        if tag.lower() == 'a' and self.current is not None:
            href, words = self.current
            self.links.append((href, ' '.join(' '.join(words).split()).lower()))
            self.current = None


def sequence_link(html, base, label):
    parser = Links()
    parser.feed(html)
    matches = []
    for href, text in parser.links:
        url = urljoin(base.rstrip('/') + '/', href)
        parts = urlparse(url)
        if (text.startswith(label.lower()) and parts.scheme == 'https'
                and parts.hostname == 'opig.stats.ox.ac.uk'
                and parts.path.endswith('.csv.gz')):
            matches.append(url)
    matches = list(dict.fromkeys(matches))
    if len(matches) != 1:
        raise ValueError('Could not uniquely resolve the sequence-only download: ' + label)
    return matches[0]


def choose_wheel(metadata):
    candidates = []
    for item in metadata.get('urls', []):
        name = item.get('filename', '')
        if (item.get('packagetype') == 'bdist_wheel' and not item.get('yanked', False)
                and '-cp313-cp313-' in name and name.endswith('.whl')
                and 'manylinux' in name and 'x86_64' in name
                and urlparse(item.get('url', '')).hostname == 'files.pythonhosted.org'):
            digest = item.get('digests', {}).get('sha256', '')
            if re.fullmatch(r'[0-9a-fA-F]{64}', digest):
                candidates.append(item)
    if not candidates:
        raise ValueError('No SHA-256-verified CPython 3.13 Linux x86-64 OpenMM wheel is listed')
    # Pinned source/version and deterministic ordering, independent of the
    # downloader computer's OS. Remote container glibc is 2.41.
    return sorted(candidates, key=lambda x: x['filename'])[0]


def failed_record(filename, role, error):
    rec = dict(filename=filename, role=role, required=True, status='failed',
               attempts=[dict(stage='source_resolution', error=str(error))])
    (ROOT / (filename + '.provenance.json')).write_text(json.dumps(rec, indent=2) + '\n', encoding='utf-8')
    message('  Could not resolve ' + filename + ': ' + str(error))
    return rec


def package(records):
    manifest = ROOT / 'manifest.json'
    manifest.write_text(json.dumps(dict(
        created_utc=utc(), fetcher_version=VERSION, files=records,
        remote_runtime='CPython 3.13; Linux x86-64; glibc 2.41',
        methodological_scope='Modern implementation of pre-2011 molecular mechanics; modern sequence datasets used only as data.',
        planned_forcefields=['Amber99SB (2006)', 'OBC generalized Born (2004)'],
        exclusions='No antibody structural-model download, machine-learning inference, protein generator, or software installation/execution.'
    ), indent=2) + '\n', encoding='utf-8')
    message('Packaging ' + str(ARCHIVE))
    with zipfile.ZipFile(ARCHIVE, 'w', zipfile.ZIP_DEFLATED, allowZip64=True) as out:
        for rec in records:
            name = rec['filename']
            if rec['status'] != 'failed':
                mode = zipfile.ZIP_STORED if name.endswith(('.gz', '.whl')) else zipfile.ZIP_DEFLATED
                out.write(ROOT / name, 'egfr_refinement_fetch/' + name, compress_type=mode)
            out.write(ROOT / (name + '.provenance.json'), 'egfr_refinement_fetch/' + name + '.provenance.json')
        out.write(manifest, 'egfr_refinement_fetch/manifest.json')
        out.write(LOG, 'egfr_refinement_fetch/download.log')
        out.write(Path(__file__).resolve(), 'fetch_egfr_refinement_inputs.py')
    return ARCHIVE


def main():
    ROOT.mkdir(parents=True, exist_ok=True)
    message('Fetching sequence-only novelty databases, a classical MM implementation, and one historical control structure.')
    message('Nothing downloaded will be installed or executed. No antibody model archives are requested.')
    records = []
    meta = fetch('openmm_pypi_metadata.json', [OPENMM_JSON], 'Pinned OpenMM release metadata and publisher SHA-256 hashes')
    records.append(meta)
    try:
        if meta['status'] == 'failed':
            raise ValueError('OpenMM release metadata download failed')
        wheel = choose_wheel(json.loads((ROOT / meta['filename']).read_text(encoding='utf-8')))
        records.append(fetch(wheel['filename'], [wheel['url']],
            'Modern implementation of classical molecular mechanics, for the remote container only',
            expected_sha=wheel['digests']['sha256'], expected_size=wheel['size']))
    except Exception as exc:
        records.append(failed_record('openmm_wheel_unresolved', 'OpenMM CPU-capable Linux wheel', exc))
    for page_name, base, items in [
        ('plabdab_landing.html', PLABDAB, [('plabdab_paired.csv.gz', 'Paired sequences'), ('plabdab_unpaired.csv.gz', 'Unpaired sequences')]),
        ('plabdab_nano_landing.html', NANO, [('plabdab_nano.csv.gz', 'All sequences')]),
    ]:
        landing = fetch(page_name, [base], 'Dataset landing page, preserving the live sequence download links')
        records.append(landing)
        for filename, label in items:
            try:
                if landing['status'] == 'failed':
                    raise ValueError('Dataset landing page download failed')
                html = (ROOT / page_name).read_text(encoding='utf-8', errors='replace')
                url = sequence_link(html, base, label)
                records.append(fetch(filename, [url], 'Antibody sequences/CDR annotations for novelty screening only; no structural models'))
            except Exception as exc:
                records.append(failed_record(filename, 'Antibody novelty sequence database', exc))
    records.append(fetch('1I1A.pdb', [
        'https://files.rcsb.org/download/1I1A.pdb',
        'https://www.ebi.ac.uk/pdbe/entry-files/download/pdb1i1a.ent'],
        '2001 FcRn-Fc acid-on binding reference control; not a binder-design seed', required=False))
    archive = package(records)
    message('DONE: ' + str(archive))
    message('Archive size: %.1f MiB' % (archive.stat().st_size / (1024 * 1024)))
    missing = [x['filename'] for x in records if x['status'] == 'failed' and x['required']]
    if missing:
        message('Important files unavailable: ' + ', '.join(missing))
    message('Upload egfr_refinement_inputs.zip even if some downloads failed; failure records are included.')
    return 0


if __name__ == '__main__':
    sys.exit(main())
