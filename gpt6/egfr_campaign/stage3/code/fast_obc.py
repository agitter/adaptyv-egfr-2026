"""Direct fixed-coordinate Amber/OBC2 interaction energies.

Classical OBC2 (2004), Debye-Huckel screening, Coulomb, Lennard-Jones,
ACE surface-area term. This is a transparent implementation, not a learned
model. Bonded and intrapartner direct terms cancel in a matched-coordinate
complex-minus-partners subtraction. No cross-partner covalent bond is allowed.

Atomic parameters are obtained from the supplied Amber99SB/OBC force field.
Units: positions/radii/sigma nm, epsilon kJ/mol, charges elementary charge.
"""
from pathlib import Path
import os,sys,json,math,time
import numpy as np
from numba import njit
R=Path('/mnt/data/egfr_campaign');S3=R/'stage3'
sys.path[:0]=[str(R/'stage2/code'),str(R/'stage2/runtime/python'),str(R/'code')]
from mm_refine import mm,app,u,FF

K=138.935456 # OpenMM electrostatic conversion; tested against the actual engine.

@njit(cache=True)
def _integral(r,ori,srj):
    if r+srj<=ori:return 0.
    U=r+srj;L=max(ori,abs(r-srj));L=max(L,1e-8)
    return .5*(1/L-1/U+.25*(r-srj*srj/r)*(1/(U*U)-1/(L*L))+.5*math.log(L/U)/r)

@njit(cache=True)
def born_pair(xyz,radius,scale,partner):
    n=len(xyz);off=radius-.009;sr=scale*off;integ=np.zeros(n);own=np.zeros(n)
    for i in range(n):
        for j in range(i):
            dx=xyz[i,0]-xyz[j,0];dy=xyz[i,1]-xyz[j,1];dz=xyz[i,2]-xyz[j,2]
            r=math.sqrt(dx*dx+dy*dy+dz*dz)
            if r<1e-7:raise ValueError('coincident atoms in OBC calculation')
            a=_integral(r,off[i],sr[j]);b=_integral(r,off[j],sr[i]);integ[i]+=a;integ[j]+=b
            if partner[i]==partner[j]:own[i]+=a;own[j]+=b
    bf=np.empty(n);bo=np.empty(n)
    for i in range(n):
        p=integ[i]*off[i];v=1/off[i]-math.tanh(p-.8*p*p+4.85*p*p*p)/radius[i];bf[i]=1/v
        p=own[i]*off[i];v=1/off[i]-math.tanh(p-.8*p*p+4.85*p*p*p)/radius[i];bo[i]=1/v
    return bf,bo

@njit(cache=True)
def interaction_coefficients(xyz,q,sigma,eps,radius,scale,partner,kappas):
    """Return LJ, SA, dielectric-independent Coulomb/GB coefficients (kJ/mol)."""
    bf,bo=born_pair(xyz,radius,scale,partner);A=0.;B=np.zeros(len(kappas));lj=0.;sa=0.;vacuum=0.;n=len(q)
    for i in range(n):
        A-=.5*K*q[i]*q[i]*(1/bf[i]-1/bo[i])
        for t in range(len(kappas)):
            kap=kappas[t];B[t]+=.5*K*q[i]*q[i]*(math.exp(-kap*bf[i])/bf[i]-math.exp(-kap*bo[i])/bo[i])
        sa+=4*math.pi*2.25936*(radius[i]+.14)**2*((radius[i]/bf[i])**6-(radius[i]/bo[i])**6)
        for j in range(i):
            dx=xyz[i,0]-xyz[j,0];dy=xyz[i,1]-xyz[j,1];dz=xyz[i,2]-xyz[j,2];r2=dx*dx+dy*dy+dz*dz;r=math.sqrt(r2)
            bb=bf[i]*bf[j];f=math.sqrt(r2+bb*math.exp(-r2/(4*bb)));c=K*q[i]*q[j]
            if partner[i]!=partner[j]:
                sr=.5*(sigma[i]+sigma[j])/r;sr6=sr**6
                lj+=4*math.sqrt(eps[i]*eps[j])*(sr6*sr6-sr6)
                vacuum+=c/r;A+=c/r-c/f
                for t in range(len(kappas)):B[t]+=c*math.exp(-kappas[t]*f)/f
            else:
                bb0=bo[i]*bo[j];f0=math.sqrt(r2+bb0*math.exp(-r2/(4*bb0)))
                A-=c*(1/f-1/f0)
                for t in range(len(kappas)):B[t]+=c*(math.exp(-kappas[t]*f)/f-math.exp(-kappas[t]*f0)/f0)
    return lj,sa,A,B,vacuum,bf,bo


def parameters(top,positions,target_chains=('A',)):
    system=FF.createSystem(top,nonbondedMethod=app.NoCutoff,constraints=None,removeCMMotion=False)
    nb=next(f for f in system.getForces() if isinstance(f,mm.NonbondedForce));gb=next(f for f in system.getForces() if isinstance(f,mm.GBSAOBCForce))
    xyz=np.asarray(positions.value_in_unit(u.nanometer),dtype=float);n=len(xyz);q=np.empty(n);sigma=np.empty(n);eps=np.empty(n);radius=np.empty(n);scale=np.empty(n)
    for i in range(n):
        a,b,c=nb.getParticleParameters(i);d,e,f=gb.getParticleParameters(i)
        q[i]=a.value_in_unit(u.elementary_charge);sigma[i]=b.value_in_unit(u.nanometer);eps[i]=c.value_in_unit(u.kilojoule_per_mole);radius[i]=e.value_in_unit(u.nanometer);scale[i]=f
        if abs(q[i]-d.value_in_unit(u.elementary_charge))>1e-8:raise ValueError('GB/direct charge mismatch')
    partner=np.array([int(a.residue.chain.id not in target_chains) for a in top.atoms()],np.int64)
    if len(set(partner))!=2:raise ValueError('Need two nonempty binding partners')
    for a,b in top.bonds():
        if partner[a.index]!=partner[b.index]:raise ValueError('Covalently linked partners cannot use this subtraction')
    return xyz,q,sigma,eps,radius,scale,partner


def evaluate(top,pos,target_chains=('A',),dielectrics=(1.,2.,4.),kappas=(0.,1.25)):
    pars=parameters(top,pos,target_chains);tic=time.time();kap=np.array(kappas,dtype=float)
    lj,sa,A,B,vac,bf,bo=interaction_coefficients(*pars,kap);states=[]
    for dielectric in dielectrics:
        for k,b in zip(kappas,B):states.append({'solute_dielectric':dielectric,'kappa_nm_inverse':k,'delta_kcal_proxy':float((lj+sa+A/dielectric+b/78.3)/4.184)})
    return {'states':states,'LJ_kcal':float(lj/4.184),'SA_kcal':float(sa/4.184),'electrostatic_A_kcal':float(A/4.184),'electrostatic_B_kcal':list(B/4.184),'vacuum_cross_kcal':float(vac/4.184),'atoms':len(pars[0]),'seconds':time.time()-tic,'interpretation':'Matched-coordinate MM/GBSA interaction proxy; excludes binding entropy and unbound conformational reorganization.'}


def engine_subtraction(top,pos,target_chains=('A',),platform='Reference'):
    def e(top,pos):
        s=FF.createSystem(top,nonbondedMethod=app.NoCutoff,constraints=None,removeCMMotion=False);integ=mm.VerletIntegrator(.001)
        ctx=mm.Context(s,integ,mm.Platform.getPlatformByName(platform));ctx.setPositions(pos);v=float(ctx.getState(getEnergy=True).getPotentialEnergy().value_in_unit(u.kilocalorie_per_mole));del ctx,integ,s;return v
    value=e(top,pos);parts=[]
    for k in [0,1]:
        m=app.Modeller(top,pos);m.delete([a for a in m.topology.atoms() if int(a.residue.chain.id not in target_chains)!=k]);parts.append(e(m.topology,m.positions))
    return value-sum(parts)


def tests():
    src=R/'stage2/intermediate/mm_binding/P00024_07_human6ARU_complex_relaxed.pdb';p=app.PDBFile(str(src));t=time.time();a=evaluate(p.topology,p.positions);b=engine_subtraction(p.topology,p.positions)
    base=a['states'][0]['delta_kcal_proxy'];err=base-b;results=[{'test':'direct OBC2 vs OpenMM Reference matched subtraction','direct':base,'openmm':b,'absolute_error_kcal':abs(err),'pass':abs(err)<.03}]
    pars=list(parameters(p.topology,p.positions));xyz=pars[0].copy();xyz[pars[-1]==1,0]+=1000.;pars[0]=xyz
    lj,sa,A,B,*_=interaction_coefficients(*pars,np.array([0.,1.25]));far=float((lj+sa+A+B[0]/78.3)/4.184)
    results.append({'test':'1000 nm separated-partner control','value_kcal':far,'pass':abs(far)<.05})
    # Partner relabeling and rigid translation invariance.
    pars=list(parameters(p.topology,p.positions));one=interaction_coefficients(*pars,np.array([0.]))[:4];pars[0]=pars[0]+[3.,-1.,5.];pars[-1]=1-pars[-1];two=interaction_coefficients(*pars,np.array([0.]))[:4]
    diff=max(float(np.max(np.abs(np.asarray(v)-np.asarray(w)))) for v,w in zip(one,two));results.append({'test':'translation and partner exchange invariance','max_difference_kJ':diff,'pass':diff<1e-6})
    record={'tests':results,'passed':all(r['pass'] for r in results),'seconds':time.time()-t,'baseline_evaluation':a}
    (S3/'reference/fast_obc_tests.json').write_text(json.dumps(record,indent=2));print(json.dumps(record,indent=2),flush=True)
    assert record['passed'],'Direct energy implementation must agree with engine before use'

if __name__=='__main__':tests()
