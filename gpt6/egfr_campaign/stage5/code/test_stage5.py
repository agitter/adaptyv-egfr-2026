"""Independent implementation checks. They do not validate biological predictions."""
from pathlib import Path
import sys,json,itertools,math
import numpy as np
S=Path('/mnt/data/egfr_campaign/stage5');R=S.parent
sys.path[:0]=[str(S/'code'),str(R/'stage4/code'),str(R/'stage3/code'),str(R/'stage2/code'),str(R/'code'),str(R/'stage2/runtime/python')]
from legacy_geometry import place,dihedral,kabsch
from build_beta_returns import extend_backwards
from continuous_carboxylate import modify
from classical_design import LIB,BB,BONDS,conformer
from build_return_loops import load_parent
from independent_mixed_priors import values,RT

def angle(a,b,c):
 x=a-b;y=c-b;return math.degrees(math.acos(float(np.clip(x@y/np.linalg.norm(x)/np.linalg.norm(y),-1,1))))
def wrapped(a,b):return abs(math.atan2(math.sin(a-b),math.cos(a-b)))
def main():
 tests=[]
 def check(name,value):tests.append({'name':name,'pass':bool(value)})
 rng=np.random.Generator(np.random.MT19937(20031990));parent,rows,tx=load_parent();right=rows[65]['atoms'];phir=dihedral(np.array(rows[64]['atoms']['C']),*[np.array(right[n]) for n in ['N','CA','C']])
 for trial in range(20):
  ph=math.radians(float(rng.uniform(-145,-110)));ps=math.radians(float(rng.uniform(115,150)));ext=extend_backwards(right,phir,ph,ps,3);chain=ext+[{n:np.array(v) for n,v in right.items()}]
  for i,r in enumerate(ext):
   nextrow=chain[i+1];check(f'beta C-N length {trial}/{i}',abs(np.linalg.norm(r['C']-nextrow['N'])-1.329)<1e-10);check(f'beta CA-C length {trial}/{i}',abs(np.linalg.norm(r['CA']-r['C'])-1.525)<1e-10);check(f'beta N-CA length {trial}/{i}',abs(np.linalg.norm(r['N']-r['CA'])-1.458)<1e-10);check(f'beta psi {trial}/{i}',wrapped(dihedral(r['N'],r['CA'],r['C'],nextrow['N']),ps)<1e-10);check(f'beta omega {trial}/{i}',wrapped(dihedral(r['CA'],r['C'],nextrow['N'],nextrow['CA']),math.pi)<1e-10);check(f'beta next phi {trial}/{i}',wrapped(dihedral(r['C'],nextrow['N'],nextrow['CA'],nextrow['C']),phir if i==2 else ph)<1e-10);check(f'beta N-CA-C angle {trial}/{i}',abs(angle(r['N'],r['CA'],r['C'])-111.2)<1e-8)
 for trial in range(20):
  aa='E';bb=np.array([rows[103]['atoms'][n] for n in BB]);a=conformer(aa,trial%len(LIB[aa]['xyz']),bb);b=modify(a,rng.uniform(-math.pi,math.pi,3));check(f'chi backbone invariant {trial}',all(np.array_equal(a[n],b[n]) for n in BB));check(f'chi bond length invariant {trial}',all(abs(np.linalg.norm(a[x]-a[y])-np.linalg.norm(b[x]-b[y]))<1e-10 for x,y in BONDS[aa]+[('CA','CB')]));check(f'chi alpha stereochemistry invariant {trial}',np.linalg.det(np.array([a[n]-a['CA'] for n in ['N','C','CB']]))*np.linalg.det(np.array([b[n]-b['CA'] for n in ['N','C','CB']]))>0)
 for trial in range(12):
  x=rng.normal(size=(18,3));q,_=np.linalg.qr(rng.normal(size=(3,3)));q[:,-1]*=np.linalg.det(q);y=x@q+rng.normal(size=3);rot,offset=kabsch(y,x);check(f'Kabsch rigid invariance {trial}',np.max(abs(y@rot+offset-x))<1e-10)
 sites=[{'kind':'histidine'},{'kind':'carboxylate','assumed_baseline_free_pKa':4.4}];states=np.array(list(itertools.product(range(3),repeat=2)));pk=[6.3,5.2];frac=[.2,.5]
 for ph in [6.5,7.4]:
  check(f'mixed constant energy pH{ph}',abs(values(states,np.full(9,-3.),sites,pk,frac,ph)[0]+3)<1e-10)
 his=np.where(states[:,0]==2,-4.,0.);acid=np.where(states[:,1]>0,3.,0.)
 check('isolated histidine acid-on',values(states,his,sites,pk,frac,7.4)[0]>values(states,his,sites,pk,frac,6.5)[0]);check('isolated unfavorable acid protonation acid-off',values(states,acid,sites,pk,frac,7.4)[0]<values(states,acid,sites,pk,frac,6.5)[0])
 for trial in range(12):
  en=rng.normal(size=9);shift=float(rng.normal());ph=float(rng.uniform(5,8));check(f'energy offset invariance {trial}',abs(values(states,en+shift,sites,pk,frac,ph)[0]-values(states,en,sites,pk,frac,ph)[0]-shift)<1e-10)
 # Read all produced sequences, including failed candidates, without claiming they all pass design gates.
 designs=[]
 for p in sorted((S/'intermediate/designs').glob('*.json')):
  d=json.loads(p.read_text())
  if 'sequence' not in d:continue
  seq=d['sequence'];cid=d['candidate_id'];check(f'standard residues/length {cid}',set(seq)<=set('ACDEFGHIKLMNPQRSTVWY') and 10<=len(seq)<=250);check(f'coordinate sequence {cid}',seq==''.join(r['aa'] for r in d['structure']));check(f'CDR slices {cid}',d['cdr_sequences']==[seq[a:b] for a,b in d['cdr_intervals_zero_based']]);check(f'no inherited ranking scores {cid}',d['metrics']['scores_valid'] is False)
  if cid!='B00000':designs.append(d)
 check('distinct stage5 candidate sequences',len({d['sequence'] for d in designs})==len(designs));check('unchanged comparator excluded from new sequence count',all(d['candidate_id']!='B00000' for d in designs))
 out={'tests':tests,'passed':sum(t['pass'] for t in tests),'failed':[t for t in tests if not t['pass']],'all_pass':all(t['pass'] for t in tests),'new_sequence_count':len(designs),'qualification':'Implementation/provenance tests only, not validation of the physical model or biological outcomes.'};(S/'reference/stage5_tests.json').write_text(json.dumps(out,indent=2));print('TESTS',out['passed'],'/',len(tests),'FAILURES',out['failed'],flush=True);assert out['all_pass']
if __name__=='__main__':main()
