"""Shrake-Rupley (1973) solvent accessibility and interface-area audit.
Bondi radii (1964), 1.4-A solvent probe, deterministic golden-angle sphere.
Uses the whole supplied receptor plus target-attached sugars for occlusion.
No binding-affinity or expression prediction is inferred from surface area.
"""
from pathlib import Path
import sys,json,math,time
import numpy as np
from scipy.spatial import cKDTree
R=Path('/mnt/data/egfr_campaign');S=R/'stage3';sys.path[:0]=[str(R/'stage2/code'),str(R/'stage2/runtime/python')]
from mm_refine import app,u
from evaluate_pool import structural_audit
RADII={'C':1.70,'N':1.55,'O':1.52,'S':1.80};PROBE=1.4

def surface(xyz,radius,partner,points=192):
    z=1.-2*(np.arange(points)+.5)/points;phi=np.arange(points)*math.pi*(3-math.sqrt(5));rho=np.sqrt(1-z*z);sphere=np.stack([rho*np.cos(phi),rho*np.sin(phi),z],axis=1);r=radius+PROBE;tree=cKDTree(xyz);a=np.zeros(len(xyz));free=a.copy()
    for i in range(len(xyz)):
        near=np.array([j for j in tree.query_ball_point(xyz[i],r[i]+r.max()) if j!=i],dtype=int);dots=xyz[i]+r[i]*sphere
        if len(near):
            blocked=(np.linalg.norm(dots[:,None,:]-xyz[near][None,:,:],axis=2)<r[near][None,:]);buried=blocked.any(1);own=blocked[:,partner[near]==partner[i]].any(1)
        else:buried=np.zeros(points,dtype=bool);own=buried.copy()
        area=4*math.pi*r[i]*r[i];a[i]=area*(1-buried.mean());free[i]=area*(1-own.mean())
    return a,free

def build(cid):
    p=app.PDBFile(str(S/'intermediate/refined'/(cid+'_human6ARU_relaxed.pdb')));xyz=np.asarray(p.positions.value_in_unit(u.angstrom));records=[];relaxed={}
    for res in p.topology.residues():
        rows=[{'xyz':xyz[a.index].tolist(),'element':a.element.symbol,'atom':a.name,'resname':res.name,'residue':int(res.id),'partner':int(res.chain.id=='B')} for a in res.atoms() if a.element!=app.element.hydrogen]
        if res.chain.id=='A':relaxed[int(res.id)+333]=rows
        else:records+=rows
    target=json.loads((R/'intermediate/human6ARU_aligned.json').read_text())
    for res in target:
        if res['human_pos'] in relaxed:records+=relaxed[res['human_pos']]
        else:
            for name,x in res['atoms'].items():records.append({'xyz':x,'element':name[0],'atom':name,'resname':res['aa'],'residue':res['human_pos'],'partner':0})
    for g in json.loads((R/'stage2/intermediate/target_glycans.json').read_text()):records.append({'xyz':g['xyz'],'element':g['atom'][0],'atom':g['atom'],'resname':g['resname'],'residue':g['resid'],'partner':0,'glycan':True})
    xyz=np.array([r['xyz'] for r in records]);radius=np.array([RADII[r['element']] for r in records]);part=np.array([r['partner'] for r in records]);return records,xyz,radius,part

def main(cid,points=192):
    start=time.time();points=int(points);records,xyz,radius,part=build(cid);complex_,free=surface(xyz,radius,part,points);burial=free-complex_;assert burial.min()>-1e-8
    rows=[]
    for pos in sorted({r['residue'] for r in records if r['partner']==1}):
        ids=[i for i,r in enumerate(records) if r['partner']==1 and r['residue']==pos];d={'position':pos,'residue':records[ids[0]]['resname'],'isolated_sasa_A2':float(sum(free[ids])),'buried_sasa_A2':float(sum(burial[ids]))};rows.append(d)
    original=json.loads((S/'intermediate/designs'/(cid+'.json')).read_text());from copy import deepcopy
    updated=deepcopy(original);brows=[r for r in p_residues(cid)];assert len(brows)==len(updated['structure'])
    for r,coords in zip(updated['structure'],brows):r['atoms']=coords
    geometry=structural_audit(updated)
    result={'candidate_id':cid,'sphere_points':points,'probe_A':PROBE,'whole_receptor_atoms_including_glycans':int(sum(part==0)),'binder_heavy_atoms':int(sum(part==1)),'buried_target_sasa_A2':float(sum(burial[part==0])),'buried_binder_sasa_A2':float(sum(burial[part==1])),'conventional_interface_area_A2':float(sum(burial)/2),'binder_total_unbound_sasa_A2':float(sum(free[part==1])),'per_binder_residue':rows,'postrefined_conformation_geometry':geometry,'seconds':time.time()-start,'qualification':'A geometric solvent-accessibility diagnostic, not a folding, affinity or expression prediction.'}
    out=S/'intermediate/surface';out.mkdir(exist_ok=True);(out/(cid+'_'+str(points)+'.json')).write_text(json.dumps(result,indent=2));print(cid,'interface_A2',result['conventional_interface_area_A2'],'alternate clashes',{n:geometry[n]['heavy_atoms_below_2A'] for n in ['1NQL','1IVO']},'seconds',round(result['seconds'],1),flush=True)

def p_residues(cid):
    p=app.PDBFile(str(S/'intermediate/refined'/(cid+'_human6ARU_relaxed.pdb')));xyz=np.asarray(p.positions.value_in_unit(u.angstrom))
    return [{a.name:xyz[a.index].tolist() for a in r.atoms() if a.element!=app.element.hydrogen and a.name!='OXT'} for r in p.topology.residues() if r.chain.id=='B']

def tests():
    xyz=np.array([[0.,0.,0.],[8.,0.,0.]]);rad=np.array([1.7,1.7]);par=np.array([0,1]);a,b=surface(xyz,rad,par,192);check=[];check.append({'test':'separated sphere has analytic area','pass':bool(np.max(abs(a-4*math.pi*3.1**2))<1e-8)});xyz[1,0]=4.;a,b=surface(xyz,rad,par,192);check.append({'test':'overlap reduces complex but not isolated areas','pass':bool(np.all(a<b) and np.all(a>=0.))});print(check,flush=True);assert all(t['pass'] for t in check);(S/'reference/surface_tests.json').write_text(json.dumps(check,indent=2))
if __name__=='__main__':
    if len(sys.argv)>1:main(*sys.argv[1:])
    else:tests()
