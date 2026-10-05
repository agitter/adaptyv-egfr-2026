"""Matched classical quench separates thermal fluctuations from persistent strain.
No extra positional/chirality force, no new trajectory, no stability inference.
"""
from pathlib import Path
import sys,json,time,hashlib
R=Path('/mnt/data/egfr_campaign');S=R/'stage7'
sys.path[:0]=[str(R/'stage4/code'),str(R/'stage2/code'),str(R/'stage2/runtime/python'),str(R/'code')]
from mm_refine import app,mm,u,FF,chirality_groups,geometry_audit
from audit_refined import audit
from Bio.SeqUtils import seq1

def main():
 start=time.time();source=Path(sys.argv[1]);tag=sys.argv[2];meta=json.loads(source.read_text());pdb=source.with_name(source.stem+'_final.pdb');out=S/'intermediate/endpoint_quench';out.mkdir(exist_ok=True)
 result=out/(tag+'_audit.json')
 if result.exists():raise FileExistsError(result)
 p=app.PDBFile(str(pdb));chains=list(p.topology.chains());assert len(chains)==1
 seq=''.join(seq1(r.name,custom_map={'HIE':'H','HID':'H','HIP':'H'}) for r in p.topology.residues());cid=meta['candidate_id']
 if cid=='NATIVE_3EAK':design={'candidate_id':cid,'sequence':seq};assert len(seq)==127
 else:design=json.loads((S/'intermediate/designs'/f'{cid}.json').read_text());assert seq==design['sequence']
 ref=chirality_groups(p.topology,p.positions);p.topology.createDisulfideBonds(p.positions)
 sysmm=FF.createSystem(p.topology,nonbondedMethod=app.NoCutoff,constraints=app.HBonds,removeCMMotion=False)
 it=mm.VerletIntegrator(.001*u.picosecond);platform=mm.Platform.getPlatformByName('CPU');platform.setPropertyDefaultValue('Threads','1')
 ctx=mm.Context(sysmm,it,platform,{'Threads':'1'});ctx.setPositions(p.positions)
 before=float(ctx.getState(getEnergy=True).getPotentialEnergy().value_in_unit(u.kilojoule_per_mole))
 mm.LocalEnergyMinimizer.minimize(ctx,10.,500)
 state=ctx.getState(getPositions=True,getEnergy=True);after=float(state.getPotentialEnergy().value_in_unit(u.kilojoule_per_mole));positions=state.getPositions()
 target=out/(tag+'_quenched.pdb')
 with target.open('w') as f:app.PDBFile.writeFile(p.topology,positions,f)
 pp=app.PDBFile(str(target));chain=next(pp.topology.chains()).id
 strict=audit(target,design,chain=chain);basic=geometry_audit(p.topology,positions,ref)
 r={'candidate_id':cid,'seed':meta['seed'],'source_had_added_chirality_restraints':meta['chirality_restraints'],'source_result':str(source),'source_result_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'source_endpoint':str(pdb),'source_endpoint_sha256':hashlib.sha256(pdb.read_bytes()).hexdigest(),'quenched_pdb':str(target),'quenched_pdb_sha256':hashlib.sha256(target.read_bytes()).hexdigest(),'forcefield':'Amber99SB2006/OBC2004','method':'500-iteration local energy minimization; HBonds constraints; no positional or added chirality restraints','initial_energy_kJ_mol':before,'final_energy_kJ_mol':after,'basic_geometry':basic,'strict_geometry':strict,'seconds':time.time()-start,'interpretation':'Phase-matched geometry diagnostic only. Energy and displacement are not folding free energies or binding affinity.'}
 result.write_text(json.dumps(r,indent=2));print(tag,'basic',basic['pass'],'strict',strict['pass'],'seconds',round(r['seconds'],1),flush=True)
if __name__=='__main__':main()
