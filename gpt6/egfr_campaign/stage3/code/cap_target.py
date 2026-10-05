"""Native-backbone peptide caps remove artificial cropped-target termini.

ACE uses the native preceding residue CA/C/O; NME uses the following residue
N/CA. No new protein fold or sequence is inferred. Binder termini are retained.
"""
from pathlib import Path
import sys,json,copy
import numpy as np
R=Path('/mnt/data/egfr_campaign');S3=R/'stage3';sys.path[:0]=[str(S3/'code'),str(R/'stage2/code'),str(R/'stage2/runtime/python')]
from mm_refine import app,u,FF,mm,geometry_audit,chirality_groups
from proton_ensemble import probe

def build(source,species,out):
    original={r['human_pos']:r for r in json.loads((R/'intermediate'/(species+'_aligned.json')).read_text())};p=app.PDBFile(str(source));xyz=np.asarray(p.positions.value_in_unit(u.angstrom));target=[r for r in p.topology.residues() if r.chain.id=='A'];binder=[r for r in p.topology.residues() if r.chain.id=='B'];assert len(target)==172
    assert 333 in original and 506 in original
    rows=[('ACE','A',0,{'CH3':original[333]['atoms']['CA'],'C':original[333]['atoms']['C'],'O':original[333]['atoms']['O']})]
    for res in target:rows.append((res.name,'A',int(res.id),{a.name:xyz[a.index] for a in res.atoms() if a.element!=app.element.hydrogen and a.name!='OXT'}))
    rows.append(('NME','A',173,{'N':original[506]['atoms']['N'],'CH3':original[506]['atoms']['CA']}))
    for res in binder:rows.append((res.name,'B',int(res.id),{a.name:xyz[a.index] for a in res.atoms() if a.element!=app.element.hydrogen}))
    lines=[];i=1;last=None
    for name,chain,rid,atoms in rows:
        if last is not None and chain!=last:lines.append('TER')
        for aname,v in atoms.items():
            x,y,z=map(float,v);lines.append(f'ATOM  {i:5d} {aname:>4s} {name:>3s} {chain}{rid:4d}    {x:8.3f}{y:8.3f}{z:8.3f}{1.:6.2f}{30.:6.2f}          {aname[0]:>2s}  ');i+=1
        last=chain
    lines+=['TER','END'];Path(out).write_text('\n'.join(lines)+'\n');q=app.PDBFile(str(out));top=q.topology;top.createDisulfideBonds(q.positions)
    capbonds=[]
    for a,b in top.bonds():
        if a.residue.name in ['ACE','NME'] or b.residue.name in ['ACE','NME']:capbonds.append([a.residue.name,a.name,b.residue.name,b.name])
    assert any(a[0]=='ACE' and a[1]=='C' and a[3]=='N' for a in capbonds) or any(a[2]=='ACE' and a[3]=='C' and a[1]=='N' for a in capbonds)
    mod=app.Modeller(top,q.positions);mod.addHydrogens(FF,pH=7.4,platform=mm.Platform.getPlatformByName('CPU'));FF.createSystem(mod.topology,nonbondedMethod=app.NoCutoff)
    audit=geometry_audit(top,q.positions,chirality_groups(top,q.positions));record={'source':str(source),'capped_pdb':str(out),'species':species,'cap_source_positions':[333,506],'cap_bonds':capbonds,'geometry_audit':audit,'forcefield_parameterization_pass':True}
    Path(str(out)+'.audit.json').write_text(json.dumps(record,indent=2));assert audit['pass'],str(audit)
    return out

def main(cid,species='human6ARU',condition='native'):
    root=S3/'intermediate/refined';stem=cid+'_'+species+('_acid' if condition=='acid' else '');source=root/(stem+'_relaxed.pdb');dest=S3/'intermediate/capped';dest.mkdir(exist_ok=True);out=dest/(stem+'_capped.pdb');prefix=dest/(stem+'_capped')
    build(source,species,out);return probe(out,prefix,('A',),{i+1:i+334 for i in range(172)})
if __name__=='__main__':main(*sys.argv[1:])
