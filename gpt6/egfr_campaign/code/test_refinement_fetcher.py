"""Offline tests: downloads are mocked. No external network or installs."""
import importlib.util,io,json,gzip,tempfile,zipfile,hashlib
from pathlib import Path
from unittest.mock import patch
spec=importlib.util.spec_from_file_location('fetcher','/mnt/data/fetch_egfr_refinement_inputs.py')
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
passed=[]
def ok(name):
 passed.append(name);print('PASS',name)
with tempfile.TemporaryDirectory() as td:
 m.ROOT=Path(td)/'fetch';m.ROOT.mkdir();m.LOG=m.ROOT/'download.log';m.ARCHIVE=Path(td)/'test.zip'
 csvdata=b'source,sequence,cdr_sequences\npaper,ACDEFGHIK,{}\n'
 gz=gzip.compress(csvdata);p=m.ROOT/'good.csv.gz';p.write_bytes(gz)
 m.validate(p,p.name);ok('valid gzip CSV')
 trunc=m.ROOT/'trunc.csv.gz';trunc.write_bytes(gz[:-6])
 try:m.validate(trunc,trunc.name);raise AssertionError('accepted truncated gzip')
 except (EOFError,OSError,ValueError):ok('truncated gzip rejected')
 bad=m.ROOT/'bad.csv.gz';bad.write_bytes(b'<html>not a dataset</html>')
 try:m.validate(bad,bad.name);raise AssertionError('accepted HTML')
 except ValueError:ok('HTML response rejected')
 html='<html><a href="static/downloads/paired.csv.gz"><b>Paired</b> sequences (.csv.gz), 11 MB</a><a href="/webapps/plabdab/static/unpaired.csv.gz">Unpaired sequences (.csv.gz)</a><a href="all.tar.gz">All PLAbDab data</a></html>'
 assert m.sequence_link(html,m.PLABDAB,'Paired sequences')==m.PLABDAB+'/static/downloads/paired.csv.gz'
 assert m.sequence_link(html,m.PLABDAB,'Unpaired sequences')=='https://opig.stats.ox.ac.uk/webapps/plabdab/static/unpaired.csv.gz'
 ok('only requested sequence links selected; relative links handled')
 try:m.sequence_link('<html><a href="https://evil.example/a.csv.gz">Paired sequences</a></html>',m.PLABDAB,'Paired sequences');raise AssertionError('accepted foreign host')
 except ValueError:ok('foreign download host rejected')
 wh=m.ROOT/'openmm-8.4.0.post2-cp313-cp313-manylinux_2_34_x86_64.whl'
 with zipfile.ZipFile(wh,'w') as z:
  z.writestr('openmm/__init__.py','# test fixture only\n')
  z.writestr('openmm-8.4.0.post2.dist-info/METADATA','Metadata-Version: 2.1\nName: OpenMM\nVersion: 8.4.0.post2\n')
 sha=m.sha256(wh);m.validate(wh,wh.name,sha,wh.stat().st_size);ok('wheel and SHA-256 validation')
 try:m.validate(wh,wh.name,'0'*64);raise AssertionError('accepted checksum mismatch')
 except ValueError:ok('publisher checksum mismatch rejected')
 metadata={'urls':[{'filename':wh.name,'packagetype':'bdist_wheel','url':'https://files.pythonhosted.org/packages/test/'+wh.name,'digests':{'sha256':sha},'size':wh.stat().st_size},
  {'filename':wh.name.replace('cp313-cp313','cp312-cp312'),'packagetype':'bdist_wheel','url':'https://files.pythonhosted.org/wrong.whl','digests':{'sha256':sha}}]}
 assert m.choose_wheel(metadata)['filename']==wh.name;ok('remote CPython ABI selected rather than local downloader ABI')
 class Response(io.BytesIO):
  def __init__(self,data,url):super().__init__(data);self.headers={'Content-Length':str(len(data))};self.url=url
  def geturl(self):return self.url
 calls=[]
 def mock_open(request,timeout=0):
  calls.append(request.full_url);return Response(gz,request.full_url)
 with patch.object(m,'urlopen',mock_open),patch.object(m.time,'sleep',lambda x:None):
  rec=m.fetch('download.csv.gz',['https://opig.stats.ox.ac.uk/test.csv.gz'],'test fixture')
  assert rec['status']=='downloaded';assert len(calls)==1
  cached=m.fetch('download.csv.gz',['https://opig.stats.ox.ac.uk/test.csv.gz'],'test fixture')
  assert cached['status']=='cached' and len(calls)==1
  ok('streaming download and cache reuse')
  badrec=m.fetch('mismatch.csv.gz',['https://opig.stats.ox.ac.uk/test.csv.gz'],'test fixture',expected_sha='0'*64)
  assert badrec['status']=='failed' and not (m.ROOT/'mismatch.csv.gz').exists()
  assert not (m.ROOT/'mismatch.csv.gz.part').exists();ok('failed download is not published or left partial')
  m.package([rec,badrec])
  with zipfile.ZipFile(m.ARCHIVE) as z:
   assert z.testzip() is None
   assert 'egfr_refinement_fetch/download.csv.gz' in z.namelist()
   assert 'egfr_refinement_fetch/mismatch.csv.gz' not in z.namelist()
   assert 'egfr_refinement_fetch/mismatch.csv.gz.provenance.json' in z.namelist()
  ok('archive preserves success and failure provenance')
print('TOTAL',len(passed),'tests passed; all network responses were fixtures')
Path('/mnt/data/egfr_campaign/reference/refinement_fetcher_test_results.json').write_text(json.dumps({'offline':True,'passed':passed,'live_network_verified':False},indent=2))
