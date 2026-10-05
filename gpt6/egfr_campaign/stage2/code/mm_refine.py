"""Amber99SB (2006), OBC (2004), L-BFGS minimization and geometry audits.
Only energy minimization is performed; no modern learned potential is used.
"""
from pathlib import Path
import sys,os,json,math,time,random,copy
import numpy as np
S=Path('/mnt/data/egfr_campaign/stage2');R=S.parent
sys.path[:0]=[str(S/'runtime/python'),str(S/'code'),str(R/'code')]
import openmm as mm
from openmm import app,unit as u
from Bio.SeqUtils import seq3
from legacy_geometry import place,dihedral
from classical_design import BONDS,BB
FF=app.ForceField('amber99sb.xml','amber99_obc.xml')

def with_terminal_oxygen(rows):
    rows=copy.deepcopy(rows)
    atoms=rows[-1]['atoms'];a={k:np.array(v) for k,v in atoms.items()}
    if 'OXT' not in atoms:
        tors=dihedral(a['N'],a['CA'],a['C'],a['O'])+math.pi
        atoms['OXT']=place(a['N'],a['CA'],a['C'],1.25,math.radians(120.8),tors).tolist()
    return rows

def write_pdb(rows,path,chain='B',terminal=True):
    if terminal:rows=with_terminal_oxygen(rows)
    lines=[];num=1
    for i,r in enumerate(rows,1):
        for name,p in r['atoms'].items():
            el=name[0];x,y,z=map(float,p)
            lines.append(f'ATOM  {num:5d} {name:>4s} {seq3(r["aa"]).upper():>3s} {chain}{i:4d}    {x:8.3f}{y:8.3f}{z:8.3f}{1.00:6.2f}{30.00:6.2f}          {el:>2s}  ');num+=1
    lines+=['TER','END'];Path(path).write_text('\n'.join(lines)+'\n')


def chirality_groups(topology,positions):
    xyz=positions.value_in_unit(u.nanometer);groups=[]
    for res in topology.residues():
        a={x.name:x.index for x in res.atoms()}
        if all(n in a for n in ['N','CA','C','CB']):groups.append(('alpha',res.index,[a[n] for n in ['N','CA','C','CB']]))
        if res.name in ['ILE','THR']:
            names=['CA','CB','CG1','CG2'] if res.name=='ILE' else ['CA','CB','OG1','CG2']
            if all(n in a for n in names):groups.append(('beta',res.index,[a[n] for n in names]))
    out=[]
    for kind,res,ids in groups:
        pts=np.array([xyz[k] for k in ids]);theta=float(dihedral(*pts));out.append((kind,res,ids,theta))
    return out


def geometry_audit(top,positions,reference_groups=None):
    from scipy.spatial import cKDTree
    coords=np.asarray(positions.value_in_unit(u.angstrom));atoms=list(top.atoms())
    graph={i:set() for i in range(len(atoms))}
    for a,b in top.bonds():graph[a.index].add(b.index);graph[b.index].add(a.index)
    heavy=[a.index for a in atoms if a.element!=app.element.hydrogen]
    excluded={}
    for i in heavy:
        seen={i};front={i}
        for _ in range(3):front={j for k in front for j in graph[k]}-seen;seen|=front
        excluded[i]=seen
    severe=[];close=[];tree=cKDTree(coords[heavy])
    for ai,bi in tree.query_pairs(2.1):
        a=heavy[ai];b=heavy[bi]
        if b in excluded[a]:continue
        d=float(np.linalg.norm(coords[a]-coords[b]));pair=[atoms[a].residue.index+1,atoms[a].name,atoms[b].residue.index+1,atoms[b].name,d]
        close.append(pair)
        if d<1.65:severe.append(pair)
    inversions=[]
    if reference_groups:
        for kind,res,ids,theta in reference_groups:
            cur=float(dihedral(*coords[ids]));delta=math.atan2(math.sin(cur-theta),math.cos(cur-theta))
            # Independent signed-volume sign check, not just restraint energy.
            if math.sin(cur)*math.sin(theta)<0:inversions.append([kind,res+1,cur,theta])
    bond_outliers=[]
    for a,b in top.bonds():
        if a.element==app.element.hydrogen or b.element==app.element.hydrogen:continue
        d=float(np.linalg.norm(coords[a.index]-coords[b.index]));upper=2.3 if (a.element==app.element.sulfur or b.element==app.element.sulfur) else 1.8
        if d<1.05 or d>upper:bond_outliers.append([a.residue.index+1,a.name,b.residue.index+1,b.name,d])
    return {'heavy_atom_count':len(heavy),'severe_nonbonded_clashes':severe,'close_nonbonded_contacts_below_2p1A':close,'stereochemical_inversions':inversions,'bond_length_outliers':bond_outliers,'pass':not(severe or inversions or bond_outliers)}


def minimize(infile,outfile,loopmask=None,max_iterations=300,obstacles=None):
    random.seed(19982003);np.random.seed(19982003)
    pdb=app.PDBFile(str(infile));mod=app.Modeller(pdb.topology,pdb.positions)
    mod.topology.createDisulfideBonds(mod.positions)
    mod.addHydrogens(FF,pH=7.4,platform=mm.Platform.getPlatformByName('CPU'))
    system=FF.createSystem(mod.topology,nonbondedMethod=app.NoCutoff,constraints=None,rigidWater=False,removeCMMotion=False)
    positions=mod.positions;ref=chirality_groups(mod.topology,positions)
    restraint=mm.CustomTorsionForce('k*(1-cos(theta-theta0))');restraint.addPerTorsionParameter('theta0');restraint.addGlobalParameter('k',1000.)
    for kind,res,ids,theta in ref:restraint.addTorsion(*ids,[theta])
    system.addForce(restraint)
    posres=mm.CustomExternalForce('kb*((x-x0)^2+(y-y0)^2+(z-z0)^2)')
    for par in ['kb','x0','y0','z0']:posres.addPerParticleParameter(par)
    nm=positions.value_in_unit(u.nanometer)
    for a in mod.topology.atoms():
        if a.name in BB:
            isloop=loopmask[a.residue.index] if loopmask is not None else False
            strength=150. if isloop else 1500.
            posres.addParticle(a.index,[strength,*nm[a.index]])
    system.addForce(posres)
    if obstacles is not None and len(obstacles):
        # Fixed glycan exclusion for binder atoms only; target glycosidic atoms
        # must not be repelled from their own covalently attached sugars.
        obs=np.array(obstacles,float)/10.
        wall=mm.CustomCompoundBondForce(1,'kg*step(r0-r)*(r0-r)^2; r=sqrt((x1-xg)^2+(y1-yg)^2+(z1-zg)^2)')
        for par in ['xg','yg','zg','r0']:wall.addPerBondParameter(par)
        wall.addGlobalParameter('kg',10000.)
        for a in mod.topology.atoms():
            if a.residue.chain.id!='B' or a.element==app.element.hydrogen:continue
            for g in obs:
                if np.linalg.norm(np.array(nm[a.index])-g)<.8:wall.addBond([a.index],[*g,.29])
        system.addForce(wall)
    for idx,f in enumerate(system.getForces()):f.setForceGroup(idx)
    integ=mm.VerletIntegrator(.001*u.picoseconds);ctx=mm.Context(system,integ,mm.Platform.getPlatformByName('CPU'),{'Threads':'1'});ctx.setPositions(positions)
    initial=float(ctx.getState(getEnergy=True).getPotentialEnergy().value_in_unit(u.kilojoule_per_mole));start=time.time()
    mm.LocalEnergyMinimizer.minimize(ctx,20.,max_iterations)
    state=ctx.getState(getEnergy=True,getPositions=True);final=float(state.getPotentialEnergy().value_in_unit(u.kilojoule_per_mole));pos=state.getPositions()
    with open(outfile,'w') as f:app.PDBFile.writeFile(mod.topology,pos,f)
    audit=geometry_audit(mod.topology,pos,ref);audit.update({'initial_energy_kjmol':initial,'final_energy_kjmol':final,'seconds':time.time()-start,'iterations_limit':max_iterations,'restraints_in_energy':True,'physical_claim':'geometry-only control; not binding or folding validation','alpha_centers':sum(x[0]=='alpha' for x in ref),'beta_centers':sum(x[0]=='beta' for x in ref)})
    # Original to relaxed backbone displacement; no superposition hides deformation.
    ids=[a.index for a in mod.topology.atoms() if a.name=='CA'];x0=np.asarray(positions.value_in_unit(u.angstrom));x1=np.asarray(pos.value_in_unit(u.angstrom));audit['ca_rms_displacement_A']=float(np.sqrt(np.mean(np.sum((x1[ids]-x0[ids])**2,axis=1))))
    del ctx,integ
    return audit

if __name__=='__main__':
    p=Path(sys.argv[1]);out=S/'intermediate/refined';out.mkdir(exist_ok=True)
    d=json.loads(p.read_text());name=d['candidate_id'];inp=out/(name+'_initial.pdb');dest=out/(name+'_relaxed.pdb');write_pdb(d['structure'],inp)
    audit=minimize(inp,dest,[r['loop'] for r in d['structure']],int(sys.argv[2]) if len(sys.argv)>2 else 300)
    audit['candidate_id']=name;(out/(name+'_audit.json')).write_text(json.dumps(audit,indent=2));print(json.dumps(audit,indent=2))
