#!/usr/bin/env python3
"""Phone-a-Friend 2: fetch sequence-only novelty data and classical MM software.

Run with Python 3 (standard library only):
    python fetch_egfr_refinement_inputs_v2.py
Upload the resulting egfr_refinement_inputs.zip.

This downloads but NEVER installs or executes any downloaded software.  The
OpenMM wheel is for the remote Linux x86-64 / CPython 3.13 computation container,
not necessarily for the computer executing this downloader.  Precomputed
antibody models and modern protein-design packages are NOT downloaded.
Version 2 adds a pinned NumPy wheel, Thera-SAbDab CSV, standalone classical
parameter files, and stricter runtime/dependency-content validation.
"""
from __future__ import annotations

import csv
from email.parser import Parser
import xml.etree.ElementTree as ET
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

VERSION = '2.0'
OPENMM_VERSION = '8.4.0.post2'
NUMPY_VERSION = '2.3.5'
ROOT = Path('egfr_refinement_fetch').resolve()
LOG = ROOT / 'download.log'
ARCHIVE = Path('egfr_refinement_inputs.zip').resolve()
MAX_DOWNLOAD = 512 * 1024 * 1024
OPENMM_JSON = 'https://pypi.org/pypi/OpenMM/' + OPENMM_VERSION + '/json'
PLABDAB = 'https://opig.stats.ox.ac.uk/webapps/plabdab'
NANO = 'https://opig.stats.ox.ac.uk/webapps/plabdab-nano'
THERA = 'https://opig.stats.ox.ac.uk/webapps/sabdab-sabpred/therasabdab/search/'
OPENMM_DATA = 'https://raw.githubusercontent.com/openmm/openmm/8.4.0/wrappers/python/openmm/app/data/'
PARAMETERS = ('amber99sb.xml', 'amber99_obc.xml', 'hydrogens.xml', 'residues.xml')
PINNED_VERSIONS = {'openmm': OPENMM_VERSION, 'numpy': NUMPY_VERSION}



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
            metadata_names = [x for x in names if x.endswith('.dist-info/METADATA')]
            if len(metadata_names) != 1:
                raise ValueError('Wheel must contain one distribution METADATA file')
            metadata = Parser().parsestr(wheel.read(metadata_names[0]).decode('utf-8'))
            project = metadata.get('Name', '').lower()
            if project not in PINNED_VERSIONS:
                raise ValueError('Unexpected wheel project: ' + project)
            if metadata.get('Version') != PINNED_VERSIONS[project]:
                raise ValueError('Wheel project version differs from the pinned version')
            if not any(x.startswith(project + '/') for x in names):
                raise ValueError('Wheel does not contain its declared package')
            if project == 'openmm':
                required_members = ['openmm/app/data/' + x for x in PARAMETERS]
                if not all(x in names for x in required_members):
                    raise ValueError('OpenMM wheel is missing required force-field/topology files')
                if not any(Path(x).name.startswith('libOpenMMCPU') and '.so' in x for x in names):
                    raise ValueError('OpenMM wheel is missing the native CPU plugin')
                if not any(Path(x).name.startswith('libOpenMM.so') for x in names):
                    raise ValueError('OpenMM wheel is missing the core native library')
            for requirement in metadata.get_all('Requires-Dist', []):
                # GPU extras are not requested. Base dependencies must be covered
                # by this offline transfer; do not silently invoke pip online.
                parts = requirement.split(';', 1)
                if len(parts) == 2 and re.search(r'\bextra\s*==', parts[1]):
                    continue
                match = re.match(r'\s*([A-Za-z0-9_.-]+)', parts[0])
                if not match or match.group(1).lower() not in PINNED_VERSIONS:
                    raise ValueError('Uncovered runtime dependency: ' + requirement)
    elif filename.endswith('.csv'):
        if b'<html' in head.lower() or b'<!doctype' in head.lower():
            raise ValueError('Received HTML instead of a CSV sequence dataset')
        with path.open(encoding='utf-8-sig', newline='') as src:
            rows = csv.reader(src)
            header = next(rows)
            if len(header) < 2 or not any('seq' in x.lower() for x in header):
                raise ValueError('CSV has no sequence column')
            if next(rows, None) is None:
                raise ValueError('CSV contains no sequence records')
    elif filename.endswith('.xml'):
        root = ET.parse(path).getroot()
        if root.tag not in ('ForceField', 'Residues'):
            raise ValueError('Unexpected force-field/topology XML root: ' + root.tag)
        if root.tag == 'ForceField' and root.find('Include') is not None:
            raise ValueError('Standalone force field has unresolved external includes')
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


def sequence_link(html, base, label, suffix='.csv.gz'):
    parser = Links()
    parser.feed(html)
    matches = []
    for href, text in parser.links:
        url = urljoin(base.rstrip('/') + '/', href)
        parts = urlparse(url)
        if (text.startswith(label.lower()) and parts.scheme == 'https'
                and parts.hostname == 'opig.stats.ox.ac.uk'
                and parts.path.endswith(suffix)):
            matches.append(url)
    matches = list(dict.fromkeys(matches))
    if len(matches) != 1:
        raise ValueError('Could not uniquely resolve the sequence-only download: ' + label)
    return matches[0]


def choose_wheel(metadata, project='openmm'):
    candidates = []
    for item in metadata.get('urls', []):
        name = item.get('filename', '')
        glibc_tags = [(int(a), int(b)) for a, b in re.findall(r'manylinux_(\d+)_(\d+)_x86_64', name)]
        compatible_glibc = any(tag <= (2, 41) for tag in glibc_tags) or 'manylinux2014_x86_64' in name
        if (item.get('packagetype') == 'bdist_wheel' and not item.get('yanked', False)
                and name.lower().startswith(project.lower() + '-') and compatible_glibc
                and '-cp313-cp313-' in name and name.endswith('.whl')
                and 'manylinux' in name and 'x86_64' in name
                and urlparse(item.get('url', '')).hostname == 'files.pythonhosted.org'):
            digest = item.get('digests', {}).get('sha256', '')
            if re.fullmatch(r'[0-9a-fA-F]{64}', digest):
                candidates.append(item)
    if not candidates:
        raise ValueError('No compatible SHA-256-verified CPython 3.13 Linux x86-64 wheel is listed for ' + project)
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
        offline_installation='OpenMM and NumPy wheels only; use --no-index --no-deps after inspecting METADATA and native libraries in the remote container.',
        runtime_checks_still_required=['Dynamic linking', 'OpenMM CPU or Reference context creation', 'Finite energies/forces', 'Amber99SB/OBC minimization on a control'],
        pH_model='Explicit fixed protonation microstates and binding polynomials, not automatic pH assignment alone; model implementation and validation remain pending.',
        coverage_limitations='Not a full public patent dump, OAS, or an exact replica of the organizers database snapshots or novelty pipeline.',
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
        out.write(Path(__file__).resolve(), Path(__file__).name)
    return ARCHIVE


def main():
    ROOT.mkdir(parents=True, exist_ok=True)
    message('Fetching offline MM runtime/dependency wheels, sequence-only novelty databases, classical parameter files, and one historical control.')
    message('Nothing downloaded will be installed or executed. No antibody model archives are requested.')
    records = []
    for project, version in [('openmm', OPENMM_VERSION), ('numpy', NUMPY_VERSION)]:
        meta = fetch(project + '_pypi_metadata.json',
            ['https://pypi.org/pypi/' + project + '/' + version + '/json'],
            'Pinned runtime metadata and publisher SHA-256 hashes')
        records.append(meta)
        try:
            if meta['status'] == 'failed':
                raise ValueError(project + ' release metadata download failed')
            wheel = choose_wheel(json.loads((ROOT / meta['filename']).read_text(encoding='utf-8')), project)
            rec = fetch(wheel['filename'], [wheel['url']],
                'Offline remote-container runtime wheel; never installed or executed by this downloader',
                expected_sha=wheel['digests']['sha256'], expected_size=wheel['size'])
            records.append(rec)
        except Exception as exc:
            records.append(failed_record(project + '_wheel_unresolved', 'Offline Linux runtime wheel', exc))
    for filename in PARAMETERS:
        records.append(fetch(filename, [OPENMM_DATA + filename],
            'Standalone classical parameter/topology data: audit and fallback; no XML scripts are executed by this downloader'))
    for page_name, base, items in [
        ('plabdab_landing.html', PLABDAB, [('plabdab_paired.csv.gz', 'Paired sequences'), ('plabdab_unpaired.csv.gz', 'Unpaired sequences')]),
        ('plabdab_nano_landing.html', NANO, [('plabdab_nano.csv.gz', 'All sequences')]),
        ('therasabdab_landing.html', THERA, [('therasabdab.csv', 'csv')]),
    ]:
        landing = fetch(page_name, [base], 'Dataset landing page, preserving the live sequence download links')
        records.append(landing)
        for filename, label in items:
            try:
                if landing['status'] == 'failed':
                    raise ValueError('Dataset landing page download failed')
                html = (ROOT / page_name).read_text(encoding='utf-8', errors='replace')
                url = sequence_link(html, base, label, '.csv.gz' if filename.endswith('.gz') else '.csv')
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
    if archive.stat().st_size > 480 * 1024 * 1024:
        message('WARNING: archive is larger than planned; retain all files and report the size before upload.')
    missing = [x['filename'] for x in records if x['status'] == 'failed' and x['required']]
    if missing:
        message('Important files unavailable: ' + ', '.join(missing))
    message('Upload egfr_refinement_inputs.zip even if some downloads failed; failure records are included.')
    return 0


if __name__ == '__main__':
    sys.exit(main())
