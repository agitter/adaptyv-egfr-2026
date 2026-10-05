from pathlib import Path
p=Path('/mnt/data/fetch_egfr_refinement_inputs.py')
s=p.read_text()
s=s.replace('python fetch_egfr_refinement_inputs.py','python fetch_egfr_refinement_inputs_v2.py')
s=s.replace("import csv\n", "import csv\nfrom email.parser import Parser\nimport xml.etree.ElementTree as ET\n")
s=s.replace("VERSION = '1.0'", "VERSION = '2.0'")
s=s.replace("OPENMM_VERSION = '8.4.0.post2'", "OPENMM_VERSION = '8.4.0.post2'\nNUMPY_VERSION = '2.3.5'")
s=s.replace("NANO = 'https://opig.stats.ox.ac.uk/webapps/plabdab-nano'", """NANO = 'https://opig.stats.ox.ac.uk/webapps/plabdab-nano'
THERA = 'https://opig.stats.ox.ac.uk/webapps/sabdab-sabpred/therasabdab/search/'
OPENMM_DATA = 'https://raw.githubusercontent.com/openmm/openmm/8.4.0/wrappers/python/openmm/app/data/'
PARAMETERS = ('amber99sb.xml', 'amber99_obc.xml', 'hydrogens.xml', 'residues.xml')
PINNED_VERSIONS = {'openmm': OPENMM_VERSION, 'numpy': NUMPY_VERSION}
""")
s=s.replace("""            if not any(x.startswith('openmm/') for x in names):
                raise ValueError('Wheel does not contain the OpenMM package')
            if not any(x.endswith('.dist-info/METADATA') for x in names):
                raise ValueError('Wheel distribution metadata is missing')""", """            metadata_names = [x for x in names if x.endswith('.dist-info/METADATA')]
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
                if len(parts) == 2 and re.search(r'\\bextra\\s*==', parts[1]):
                    continue
                match = re.match(r'\\s*([A-Za-z0-9_.-]+)', parts[0])
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
            raise ValueError('Standalone force field has unresolved external includes')""")
# More conservative wheel filtering and reusable selector.
s=s.replace('def choose_wheel(metadata):', "def choose_wheel(metadata, project='openmm'):")
s=s.replace("""        if (item.get('packagetype') == 'bdist_wheel' and not item.get('yanked', False)
                and '-cp313-cp313-' in name and name.endswith('.whl')""", """        glibc_tags = [(int(a), int(b)) for a, b in re.findall(r'manylinux_(\\d+)_(\\d+)_x86_64', name)]
        compatible_glibc = any(tag <= (2, 41) for tag in glibc_tags) or 'manylinux2014_x86_64' in name
        if (item.get('packagetype') == 'bdist_wheel' and not item.get('yanked', False)
                and name.lower().startswith(project.lower() + '-') and compatible_glibc
                and '-cp313-cp313-' in name and name.endswith('.whl')""")
s=s.replace("raise ValueError('No SHA-256-verified CPython 3.13 Linux x86-64 OpenMM wheel is listed')", "raise ValueError('No compatible SHA-256-verified CPython 3.13 Linux x86-64 wheel is listed for ' + project)")
# Generalized sequence URL selection for Thera CSV only (not spreadsheets).
s=s.replace("def sequence_link(html, base, label):", "def sequence_link(html, base, label, suffix='.csv.gz'):")
s=s.replace("and parts.path.endswith('.csv.gz')):", "and parts.path.endswith(suffix)):")
# Preserve actual dependency declarations as data, but do not install anything.
s=s.replace("""        exclusions='No antibody structural-model download, machine-learning inference, protein generator, or software installation/execution.'""", """        offline_installation='OpenMM and NumPy wheels only; use --no-index --no-deps after inspecting METADATA and native libraries in the remote container.',
        runtime_checks_still_required=['Dynamic linking', 'OpenMM CPU or Reference context creation', 'Finite energies/forces', 'Amber99SB/OBC minimization on a control'],
        pH_model='Explicit fixed protonation microstates and binding polynomials, not automatic pH assignment alone; model implementation and validation remain pending.',
        coverage_limitations='Not a full public patent dump, OAS, or an exact replica of the organizers database snapshots or novelty pipeline.',
        exclusions='No antibody structural-model download, machine-learning inference, protein generator, or software installation/execution.'""")
s=s.replace("out.write(Path(__file__).resolve(), 'fetch_egfr_refinement_inputs.py')", "out.write(Path(__file__).resolve(), Path(__file__).name)")
s=s.replace("message('Fetching sequence-only novelty databases, a classical MM implementation, and one historical control structure.')", "message('Fetching offline MM runtime/dependency wheels, sequence-only novelty databases, classical parameter files, and one historical control.')")
start=s.index("    meta = fetch('openmm_pypi_metadata.json'")
end=s.index("    for page_name, base, items in [",start)
s=s[:start]+"""    for project, version in [('openmm', OPENMM_VERSION), ('numpy', NUMPY_VERSION)]:
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
"""+s[end:]
s=s.replace("""        ('plabdab_nano_landing.html', NANO, [('plabdab_nano.csv.gz', 'All sequences')]),""", """        ('plabdab_nano_landing.html', NANO, [('plabdab_nano.csv.gz', 'All sequences')]),
        ('therasabdab_landing.html', THERA, [('therasabdab.csv', 'csv')]),""")
s=s.replace("url = sequence_link(html, base, label)", "url = sequence_link(html, base, label, '.csv.gz' if filename.endswith('.gz') else '.csv')")
s=s.replace("    message('Archive size: %.1f MiB' % (archive.stat().st_size / (1024 * 1024)))", """    message('Archive size: %.1f MiB' % (archive.stat().st_size / (1024 * 1024)))
    if archive.stat().st_size > 480 * 1024 * 1024:
        message('WARNING: archive is larger than planned; retain all files and report the size before upload.')""")
# Replace obsolete no-model text in docstring with explicit additions.
s=s.replace('antibody models and modern protein-design packages are NOT downloaded.', 'antibody models and modern protein-design packages are NOT downloaded.\nVersion 2 adds a pinned NumPy wheel, Thera-SAbDab CSV, standalone classical\nparameter files, and stricter runtime/dependency-content validation.')
Path('/mnt/data/fetch_egfr_refinement_inputs_v2.py').write_text(s)
Path('/mnt/data/egfr_campaign/code/fetch_egfr_refinement_inputs_v2.py').write_text(s)
print('Created v2:',len(s),'bytes')
