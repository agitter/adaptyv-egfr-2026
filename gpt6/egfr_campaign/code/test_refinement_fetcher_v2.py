"""Offline downloader tests; fabricated HTTP responses, no scientific engine run."""
import importlib.util,io,json,gzip,tempfile,zipfile,hashlib,unittest
from pathlib import Path
from unittest.mock import patch
spec=importlib.util.spec_from_file_location('fetcher','/mnt/data/fetch_egfr_refinement_inputs_v2.py')
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)

def wheel_bytes(project='openmm',omit=None,dependencies=None,version=None):
    buf=io.BytesIO();version=version or m.PINNED_VERSIONS[project]
    members={project+'/__init__.py':'# fixture only; never executed\n'}
    meta='Metadata-Version: 2.1\nName: '+project+'\nVersion: '+version+'\n'
    deps=dependencies if dependencies is not None else (['numpy'] if project=='openmm' else [])
    meta+=''.join('Requires-Dist: '+x+'\n' for x in deps)
    members[project+'-'+version+'.dist-info/METADATA']=meta
    if project=='openmm':
        members.update({'openmm/app/data/'+x:'<ForceField></ForceField>\n' for x in m.PARAMETERS})
        members['OpenMM.libs/lib/libOpenMM.so']='fixture native library, not runnable\n'
        members['OpenMM.libs/lib/plugins/libOpenMMCPU.so']='fixture native library, not runnable\n'
    if omit: members.pop(omit,None)
    with zipfile.ZipFile(buf,'w') as z:
        for name,data in members.items(): z.writestr(name,data)
    return buf.getvalue()

def file_item(project='openmm',abi='cp313',platform='manylinux_2_34_x86_64',**kw):
    data=wheel_bytes(project)
    name=project+'-'+m.PINNED_VERSIONS[project]+'-'+abi+'-'+abi+'-'+platform+'.whl'
    out={'filename':name,'packagetype':'bdist_wheel','url':'https://files.pythonhosted.org/packages/test/'+name,
         'digests':{'sha256':hashlib.sha256(data).hexdigest()},'size':len(data)}
    out.update(kw);return out

class Response(io.BytesIO):
    def __init__(self,data,url):
        super().__init__(data);self.url=url;self.headers={'Content-Length':str(len(data))}
    def geturl(self):return self.url

class Tests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup)
        m.ROOT=Path(self.temp.name)/'fetch';m.ROOT.mkdir();m.LOG=m.ROOT/'download.log';m.ARCHIVE=Path(self.temp.name)/'out.zip'
        self.csv=b'source,sequence,cdr_sequences\npaper,ACDEFGHIK,{}\n';self.gz=gzip.compress(self.csv)
    def validate(self,name,data,**kw):
        p=m.ROOT/name;p.write_bytes(data);return m.validate(p,name,**kw)
    def test_01_gzip(self):self.validate('data.csv.gz',self.gz)
    def test_02_truncated_gzip(self):
        with self.assertRaises((EOFError,OSError,ValueError)):self.validate('data.csv.gz',self.gz[:-6])
    def test_03_html_instead_of_gzip(self):
        with self.assertRaises(ValueError):self.validate('data.csv.gz',b'<html>failure message</html>')
    def test_04_csv(self):self.validate('thera.csv',self.csv)
    def test_05_html_instead_of_csv(self):
        with self.assertRaises(ValueError):self.validate('thera.csv',b'<html>failure message</html>')
    def test_06_no_sequence_csv(self):
        with self.assertRaises(ValueError):self.validate('thera.csv',b'aaa,bbbb,cccc,dddd\n1,2,3,4\n')
    def test_07_empty_csv(self):
        with self.assertRaises(ValueError):self.validate('thera.csv',b'source,sequence,cdr_sequences\n')
    def test_08_forcefield(self):self.validate('amber.xml',b'<ForceField><Residues/></ForceField>\n')
    def test_09_hydrogen_definitions(self):self.validate('hydrogens.xml',b'<Residues><Residue name="ALA"/></Residues>\n')
    def test_10_wrong_xml(self):
        with self.assertRaises(ValueError):self.validate('amber.xml',b'<html><body>failure</body></html>\n')
    def test_11_unresolved_include(self):
        with self.assertRaises(ValueError):self.validate('amber.xml',b'<ForceField><Include file="missing.xml"/></ForceField>\n')
    def test_12_openmm_wheel(self):self.validate('openmm.whl',wheel_bytes())
    def test_13_numpy_wheel(self):self.validate('numpy.whl',wheel_bytes('numpy'))
    def test_14_wheel_no_cpu(self):
        with self.assertRaises(ValueError):self.validate('openmm.whl',wheel_bytes(omit='OpenMM.libs/lib/plugins/libOpenMMCPU.so'))
    def test_15_wheel_no_forcefield(self):
        with self.assertRaises(ValueError):self.validate('openmm.whl',wheel_bytes(omit='openmm/app/data/amber99sb.xml'))
    def test_16_dependency_missing(self):
        with self.assertRaises(ValueError):self.validate('openmm.whl',wheel_bytes(dependencies=['numpy','surprise_dependency']))
    def test_17_gpu_extra_not_required(self):self.validate('openmm.whl',wheel_bytes(dependencies=['numpy','nvidia-driver; extra == "cuda12"']))
    def test_18_wrong_wheel_version(self):
        with self.assertRaises(ValueError):self.validate('openmm.whl',wheel_bytes(version='8.5.0'))
    def test_19_wheel_checksum(self):
        with self.assertRaises(ValueError):self.validate('openmm.whl',wheel_bytes(),expected_sha='0'*64)
    def test_20_wheel_abi(self):
        good=file_item();bad=file_item(abi='cp312')
        self.assertEqual(m.choose_wheel({'urls':[bad,good]})['filename'],good['filename'])
    def test_21_wheel_glibc(self):
        with self.assertRaises(ValueError):m.choose_wheel({'urls':[file_item(platform='manylinux_2_99_x86_64')]})
    def test_22_wheel_arm(self):
        with self.assertRaises(ValueError):m.choose_wheel({'urls':[file_item(platform='manylinux_2_34_aarch64')]})
    def test_23_numpy_selector(self):
        item=file_item('numpy');self.assertEqual(m.choose_wheel({'urls':[item]},'numpy')['filename'],item['filename'])
    def test_24_sequence_links(self):
        html='<html><a href="static/paired.csv.gz"><b>Paired</b> sequences (csv.gz)</a><a href="models.tar.gz">Models</a></html>'
        self.assertEqual(m.sequence_link(html,m.PLABDAB,'Paired sequences'),m.PLABDAB+'/static/paired.csv.gz')
    def test_25_thera_csv_link(self):
        html='<html><a href="/static/thera.xlsx">xlsx</a> / <a href="/static/thera.csv">csv</a></html>'
        self.assertEqual(m.sequence_link(html,m.THERA,'csv','.csv'),'https://opig.stats.ox.ac.uk/static/thera.csv')
    def test_26_foreign_link(self):
        with self.assertRaises(ValueError):m.sequence_link('<a href="https://evil.example/a.csv">csv</a>',m.THERA,'csv','.csv')
    def test_27_download_cache(self):
        calls=[]
        def mock(req,timeout):calls.append(req.full_url);return Response(self.gz,req.full_url)
        with patch.object(m,'urlopen',mock):
            a=m.fetch('data.csv.gz',['https://example.org/data.csv.gz'],'test')
            b=m.fetch('data.csv.gz',['https://example.org/data.csv.gz'],'test')
        self.assertEqual(a['status'],'downloaded');self.assertEqual(b['status'],'cached');self.assertEqual(len(calls),1)
    def test_28_bad_download_cleanup(self):
        with patch.object(m,'urlopen',lambda req,timeout:Response(self.gz,req.full_url)),patch.object(m.time,'sleep',lambda x:None):
            rec=m.fetch('bad.csv.gz',['https://example.org/bad.csv.gz'],'test',expected_sha='0'*64)
        self.assertEqual(rec['status'],'failed');self.assertFalse((m.ROOT/'bad.csv.gz.part').exists())
    def test_29_full_workflow_and_cache(self):
        urls={}
        for project in m.PINNED_VERSIONS:
            item=file_item(project);urls[item['url']]=wheel_bytes(project)
            urls['https://pypi.org/pypi/'+project+'/'+m.PINNED_VERSIONS[project]+'/json']=json.dumps({'urls':[item]}).encode()
        for name in m.PARAMETERS:urls[m.OPENMM_DATA+name]=b'<ForceField><Residues/></ForceField>\n'
        for base,items in [(m.PLABDAB,[('paired.csv.gz','Paired sequences'),('unpaired.csv.gz','Unpaired sequences')]),
                           (m.NANO,[('nano.csv.gz','All sequences')]),(m.THERA,[('thera.csv','csv')])]:
            anchors=[]
            for filename,label in items:
                url=base.rstrip('/')+'/static/'+filename
                urls[url]=self.gz if filename.endswith('.gz') else self.csv
                anchors.append('<a href="'+url+'">'+label+'</a>')
            urls[base]=('<html>'+''.join(anchors)+'</html>').encode()
        urls['https://files.rcsb.org/download/1I1A.pdb']=b'HEADER    HISTORICAL CONTROL FIXTURE NOT REAL PDB\n'
        calls=[]
        def mock(req,timeout):calls.append(req.full_url);return Response(urls[req.full_url],req.full_url)
        with patch.object(m,'urlopen',mock),patch.object(m.time,'sleep',lambda x:None):
            self.assertEqual(m.main(),0);n=len(calls);self.assertEqual(m.main(),0);self.assertEqual(len(calls),n)
        manifest=json.loads((m.ROOT/'manifest.json').read_text())
        self.assertEqual(len(manifest['files']),16);self.assertTrue(all(r['status']=='cached' for r in manifest['files']))
        with zipfile.ZipFile(m.ARCHIVE) as z:
            self.assertIsNone(z.testzip());self.assertIn('egfr_refinement_fetch/therasabdab.csv',z.namelist())
            self.assertIn('fetch_egfr_refinement_inputs_v2.py',z.namelist())

if __name__=='__main__':
    suite=unittest.defaultTestLoader.loadTestsFromTestCase(Tests)
    result=unittest.TextTestRunner(verbosity=2).run(suite)
    out={'tests_run':result.testsRun,'failures':len(result.failures),'errors':len(result.errors),'offline':True,
         'live_downloads_verified':False,'openmm_runtime_tested':False,'scientific_results_generated':False}
    Path('/mnt/data/egfr_campaign/reference/refinement_fetcher_v2_test_results.json').write_text(json.dumps(out,indent=2)+'\n')
    raise SystemExit(not result.wasSuccessful())
