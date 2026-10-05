import sys,json,math,numpy as np
from classical_design import *

def run():
    results=[]
    def check(name,ok):
        assert bool(ok),name
        results.append(name)
    check('unchanged proton affinity gives no pH effect',np.max(np.abs(proton_free_energy(np.zeros(5),6.5)))<1e-12)
    check('stabilized protonation gives acid-on',proton_free_energy(-3.,6.5)<proton_free_energy(-3.,7.4))
    check('destabilized protonation gives acid-off',proton_free_energy(3.,6.5)>proton_free_energy(3.,7.4))
    for pka in [5.5,6.,6.5,7.,7.5]:
        for shift in [-12.,-5.,-1.,0.,1.,5.]:
            a=float(proton_free_energy(shift,6.5,pka));b=float(proton_free_energy(shift,7.4,pka));check(f'linkage bound {pka} {shift}',abs(b-a)<=RT*math.log(10)*.9+1e-12)
            check(f'compiled equivalence {pka} {shift}',abs(a-proton_energy_numba(shift,6.5,pka))<1e-12)
    for f in ['human6ARU_chainA.json','3EAK_chainA.json']:
        rows=json.loads((O/f).read_text())
        signs=[]
        for r in rows:
            a=r['atoms']
            if all(k in a for k in ['N','CA','C','CB']):signs.append(np.linalg.det(np.array([np.array(a[k])-a['CA'] for k in ['N','C','CB']])))
        check('native L chirality '+f,all(q>0 for q in signs))
    bb=np.array([[0.,0.,0.],[1.45,0,0],[2.,1.42,0],[2.,2.5,0]])
    for aa in AA:
        for k in range(len(LIB[aa]['xyz'])):
            atoms=conformer(aa,k,bb)
            if aa!='G':check(f'constructed chirality {aa} {k}',np.linalg.det(np.array([atoms[x]-atoms['CA'] for x in ['N','C','CB']]))>0)
    # Independent scalar brute-force hydrogen-free pair test and symmetry.
    atoms=conformer('S',0,bb);a=atom_data('S',atoms,['CB','OG']);b=a.copy();b[:,:3]+=np.array([8.,0,0])
    check('pair symmetry',np.allclose(pair_score(a,b),pair_score(b,a)))
    far=a.copy();far[:,:3]+=50
    check('no long-range vdW/hbond artifact',np.allclose(pair_score(a,far),0))
    check('overlap strongly penalized',pair_score(a,a)[0]>100)
    check('opposite charges attractive',screened_charge(-1,1,4)<0)
    check('like charges repel',screened_charge(1,1,4)>0)
    (S/'reference/classical_tests.json').write_text(json.dumps({'passed':len(results),'tests':results},indent=2))
    print('PASSED',len(results))
if __name__=='__main__':run()
