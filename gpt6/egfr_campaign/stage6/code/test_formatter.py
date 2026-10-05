"""Offline input-validation tests for the optional submission formatter."""
from pathlib import Path
import subprocess,sys,csv,json,hashlib
S=Path('/mnt/data/egfr_campaign/stage6');T=S/'intermediate/formatter_tests';T.mkdir(exist_ok=True)
script=S/'code/make_submission.py';a='ACDEFGHIKLMNPQRSTVWY'
rows=[(f'TEST_{i:03}', 'ACDEFGHIKL'+a[i//20]+a[i%20]) for i in range(100)]
fasta=T/'synthetic_FORMAT_TEST_ONLY.fasta';fasta.write_text(''.join('>'+n+'\n'+s+'\n' for n,s in rows));results=[]
def run(label,args,code):
 p=subprocess.run([sys.executable,str(script),*map(str,args)],text=True,capture_output=True);r={'test':label,'expected_returncode':code,'actual_returncode':p.returncode,'stdout':p.stdout,'stderr':p.stderr,'pass':p.returncode==code};results.append(r);assert r['pass'],r
 return p
out=T/'baseline.csv';run('baseline20',['--fasta',fasta,'--output',out],0)
with out.open() as f:r=list(csv.DictReader(f))
assert [(x['name'],x['sequence']) for x in r]==rows[:20]
reject=T/'reject.txt';reject.write_text('# Only test records\nTEST_000\nTEST_001\nTEST_002\n')
run('reserve_replacement',['--fasta',fasta,'--reject-file',reject,'--output',T/'revised.csv'],0)
with (T/'revised.csv').open() as f:r=list(csv.DictReader(f))
assert [(x['name'],x['sequence']) for x in r]==rows[3:23]
run('track3_cap',['--fasta',fasta,'--count','21','--output',T/'should_not_exist.csv'],2)
reject.write_text('UNKNOWN_IDENTIFIER\n');run('unknown_rejection',['--fasta',fasta,'--reject-file',reject,'--output',T/'bad.csv'],2)
reject.write_text('\n'.join(n for n,s in rows));run('insufficient_reserves',['--fasta',fasta,'--reject-file',reject,'--output',T/'bad.csv'],2)
for label,text in [('invalid_amino_acid','>BAD\nACDEFGHIKLX\n'),('short_sequence','>BAD\nACD\n'),('duplicate_sequence','>A\nACDEFGHIKLM\n>B\nACDEFGHIKLM\n'),('duplicate_identifier','>A\nACDEFGHIKLM\n>A\nACDEFGHIKLN\n'),('invalid_identifier','>A instruction\nACDEFGHIKLM\n'),('header_missing','ACDEFGHIKLM\n'),('empty','')]:
 p=T/(label+'.fasta');p.write_text(text);run(label,['--fasta',p,'--output',T/'bad.csv'],2)
h=hashlib.sha256(fasta.read_bytes()).hexdigest();run('prevent_input_overwrite',['--fasta',fasta,'--output',fasta],2);assert hashlib.sha256(fasta.read_bytes()).hexdigest()==h
(S/'reference/formatter_tests.json').write_text(json.dumps({'passed':len(results),'count':len(results),'results':results,'synthetic_inputs':'FORMAT TESTS ONLY, not candidate protein designs.'},indent=2));print('FORMATTER TESTS',len(results),'PASSED')
