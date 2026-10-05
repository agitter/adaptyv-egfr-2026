"""Independent final-file, source-coordinate, and process-state checks."""
from pathlib import Path
import collections,csv,hashlib,json,os,shlex,subprocess,sys
from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1
R=Path('/mnt/data/egfr_campaign');S=R/'stage6';O=S/'output'
results=[]
def check(name,value,details=None):
 r={'check':name,'pass':bool(value)}
 if details is not None:r['details']=details
 results.append(r)
 if not value:raise AssertionError(r)
def fasta(path):
 rows=[];n=None;ss=[]
 for line in Path(path).read_text().splitlines():
  if line.startswith('>'):
   if n is not None:rows.append((n,''.join(ss)))
   n=line[1:];ss=[]
  elif line:ss.append(line)
 if n is not None:rows.append((n,''.join(ss)))
 return rows
panel=json.loads((O/'final_panel_private.json').read_text());fa=fasta(O/'egfr_ranked_100.fasta')
check('FASTA_100_unique',len(fa)==100 and len({n for n,s in fa})==100 and len({s for n,s in fa})==100)
check('FASTA_matches_private_records',fa==[(r['public_name'],r['sequence']) for r in panel])
check('format_standard_AA_and_length',all(10<=len(s)<=250 and set(s)<=set('ACDEFGHIKLMNPQRSTVWY') for n,s in fa))
with (O/'egfr_track3_top20.csv').open(newline='') as f:
 reader=csv.DictReader(f);columns=reader.fieldnames;rows=list(reader)
check('CSV_required_columns_only',columns==['name','sequence','molecule_class'])
check('CSV_top20_order_and_sequence',len(rows)==20 and [(r['name'],r['sequence']) for r in rows]==fa[:20])
check('CSV_all_nanobody',all(r['molecule_class']=='nanobody' for r in rows))
check('no_functional_success_labels',all(set(r['biological_criteria'].values())=={'UNVERIFIED'} for r in panel))
check('all_official_novelty_pending',all(r['official_novelty_status']=='PENDING_ADAPTYV' for r in panel))
check('local_CDR3_gate',all(r['local_novelty']['max_edit_identity']<.7 for r in panel))
check('top20_two_species_detailed_geometry',all(r['human_models'] and r['mouse_models'] for r in panel[:20]))
check('glycan_results_explicit',all(r['glycan_stress'] and r['glycan_stress']['conservative_panel_filter_pass'] for r in panel))
for r in panel:
 p=Path(r['design_source']);d=json.loads(p.read_text())
 check('source_'+r['public_name'],hashlib.sha256(p.read_bytes()).hexdigest()==r['design_sha256'] and d['sequence']==r['sequence'] and ''.join(v['aa'] for v in d['structure'])==r['sequence'])
# Read actual complex coordinate files independently of the modeling library.
parser=PDBParser(QUIET=True);seen={};nmodels=0
maps={'HID':'H','HIE':'H','HIP':'H','HIN':'H','CYX':'C','CYM':'C','ASH':'D','GLH':'E','LYN':'K'}
for r in panel:
 for m in r['human_models']+r['mouse_models']:
  p=Path(m['coordinate_source']);key=(str(p),r['sequence'])
  if key in seen:continue
  model=parser.get_structure('check',str(p))[0]
  check('complex_chain_B_present_'+p.stem,'B' in model)
  seq=''.join(seq1(x.resname,custom_map=maps) for x in model['B'].get_residues())
  check('coordinate_sequence_'+p.stem,seq==r['sequence'],{'coordinate_file':str(p),'sequence_length':len(seq)})
  seen[key]=True;nmodels+=1
check('all_top20_models_are_positive_human_direction',all(m['independent_priors']['all_acid_on'] for r in panel[:20] for m in r['human_models']))
# The correction utility must reproduce the original CSV with no exclusions.
t=S/'intermediate/formatter_tests';t.mkdir(exist_ok=True)
q=t/'real_panel_roundtrip.csv';p=subprocess.run([sys.executable,str(O/'make_submission.py'),'--fasta',str(O/'egfr_ranked_100.fasta'),'--output',str(q)],text=True,capture_output=True)
check('real_panel_formatter_exit',p.returncode==0,{'stdout':p.stdout,'stderr':p.stderr})
check('real_panel_formatter_byte_identical',q.read_bytes()==(O/'egfr_track3_top20.csv').read_bytes())
# Verify no outstanding scientific worker or BLAST process is left behind.
raw=subprocess.check_output(['ps','-eo','pid,args'],text=True);running=[]
for line in raw.splitlines()[1:]:
 try:
  pidtext,cmd=line.strip().split(None,1);pid=int(pidtext);args=shlex.split(cmd)
 except (ValueError,IndexError):continue
 if pid==os.getpid() or not args:continue
 executable=Path(args[0]).name
 if executable.startswith('python') and len(args)>1 and args[1].startswith(str(S/'code')):running.append({'pid':pid,'command':args})
 elif args[0].startswith(str(S/'runtime/bin')):running.append({'pid':pid,'command':args})
check('no_remaining_campaign_workers',not running,running)
report={'checks':len(results),'passed':sum(r['pass'] for r in results),'coordinate_models_independently_sequence_checked':nmodels,'active_workers':running,'results':results}
(S/'reference/final_postflight.json').write_text(json.dumps(report,indent=2));print('POSTFLIGHT',len(results),'checks passed;',nmodels,'coordinate sequences verified')
