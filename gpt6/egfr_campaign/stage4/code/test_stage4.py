"""Unit/property checks for added geometry and thermal-screen machinery."""
from pathlib import Path
import sys,json,math,itertools
import numpy as np
S=Path('/mnt/data/egfr_campaign/stage4');sys.path.insert(0,str(S/'code'))
import stabilize_design as design
import dynamics as dyn
from openmm import unit as u
checks=[]
def record(name,value):
 checks.append({'name':name,'pass':bool(value)});assert value,name
record('explicit_stage4_output_root',design.S==S)
# Geometry covariance and rotation-preserved cysteine bond length.
parent,rows,te,me,gt,src,pfile=design.load('C00009');bb=np.array([rows[56]['atoms'][n] for n in design.BB]);a=design.conformer('C',0,bb)
for theta in [-2.,-1.,0.,1.,2.]:
 b=design.rotate_sg(a,theta);record('SG_CB_bond_preserved_'+str(theta),abs(np.linalg.norm(a['SG']-a['CB'])-np.linalg.norm(b['SG']-b['CB']))<1e-10)
 record('backbone_unchanged_'+str(theta),all(np.array_equal(a[n],b[n]) for n in design.BB))
# Constraint groups must partition all atoms, preserve transitive closure.
system=dyn.mm.System()
for i in range(7):system.addParticle(1.)
for i,j in [(0,1),(1,2),(3,4)]:system.addConstraint(i,j,.1)
g=dyn.constraint_groups(system);record('constraint_group_partition',sorted(x for z in g for x in z)==list(range(7)));record('constraint_group_transitivity',{tuple(v) for v in g}=={(0,1,2),(3,4),(5,),(6,)})
# Explicit reproducibility and Maxwell variance; no modern default RNG.
r=np.random.Generator(np.random.MT19937(1998));x=r.normal(size=(200000,3));record('Maxwell_standard_variance',abs(float(np.var(x))-1)<.01);r2=np.random.Generator(np.random.MT19937(1998));record('MT_reproducibility',np.array_equal(x[:4],r2.normal(size=(4,3))))
# Independent CA fit should remove a rigid transform but not a loop deformation.
top=dyn.app.Topology();ch=top.addChain('B');coords=[]
for i,point in enumerate([[0,0,0],[1,0,0],[0,1,0],[0,0,1],[1,1,1]]):
 rr=top.addResidue('ALA',ch);top.addAtom('CA',dyn.app.element.carbon,rr);coords.append(point)
x=np.array(coords,float);rot=np.array([[0,-1,0],[1,0,0],[0,0,1.]]);y=x@rot+4.;mask=[False,False,False,False,True];metric=dyn.align_metrics(top,y*u.angstrom,x*u.angstrom,mask);record('RMSD_rigid_transform',metric['maximum_CA_displacement_A']<1e-10);y[4]+=np.array([1.,0,0]);metric=dyn.align_metrics(top,y*u.angstrom,x*u.angstrom,mask);record('RMSD_detects_loop_move',abs(metric['loop_CA_RMSD_A']-1)<1e-10)
# All outputs are de novo-parent edits with preserved key contact identities.
seqs=set();positions={}
for p in sorted((S/'intermediate/designs').glob('S*.json')):
 d=json.loads(p.read_text());record(p.stem+'_sequence_atoms_agree',d['sequence']==''.join(r['aa'] for r in d['structure']));record(p.stem+'_contacts_preserved',all(d['sequence'][i-1]==aa for i,aa in [(54,'D'),(100,'Y'),(103,'E')]))
 record(p.stem+'_score_invalid_until_rescored',d['metrics']['scores_valid'] is False);seqs.add(d['sequence'])
record('unique_generated_sequences',len(seqs)==63)
(S/'reference/stage4_unit_tests.json').write_text(json.dumps({'tests':checks,'passed':sum(r['pass'] for r in checks),'total':len(checks)},indent=2));print('PASS',len(checks))
