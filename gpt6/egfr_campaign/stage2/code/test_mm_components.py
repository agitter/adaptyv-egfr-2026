"""Independent zero-interaction limit test for matched-component MM/GBSA."""
from mm_binding_probe import *

def main():
    p=S/'intermediate/mm_binding/P00024_07_human6ARU_complex_relaxed.pdb'
    pdb=app.PDBFile(str(p));mod=app.Modeller(pdb.topology,pdb.positions)
    mod.delete([a for a in mod.topology.atoms() if a.element==app.element.hydrogen])
    x=np.array(mod.positions.value_in_unit(u.nanometer))
    for a in mod.topology.atoms():
        if a.residue.chain.id=='B':x[a.index]+=np.array([1000.,0.,0.])
    mod.positions=x*u.nanometer
    variants=['HIE' if r.name=='HIS' else None for r in mod.topology.residues()]
    e=mmgbsa(mod,variants)
    result={'test':'matched components approach zero interaction at 1000 nm separation',
      'separation_nm':1000.,'delta_kcal_proxy':e['delta_kcal_proxy'],
      'tolerance_kcal':.20,'pass':abs(e['delta_kcal_proxy'])<.20,
      'interpretation':'Tests subtraction consistency only; does not validate affinity prediction.'}
    (S/'reference/mm_component_test.json').write_text(json.dumps(result,indent=2))
    print(json.dumps(result,indent=2))
    assert result['pass'],result
if __name__=='__main__':main()
