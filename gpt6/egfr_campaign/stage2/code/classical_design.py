"""Classical discrete side-chain design. No learned models or sequence seeds.

Geometric model: empirical conformers, pairwise softened van der Waals and
angular hydrogen bonding, screened Coulomb formal-charge interactions.
Binding-polynomial pH averaging is applied to protonatable histidines.
All coefficients are declared search heuristics, NOT calibrated binding energies.
Sequence search uses Metropolis simulated annealing and MT19937.
"""
from pathlib import Path
import json,sys,math,copy
import numpy as np
from scipy.spatial.distance import cdist
from scipy.special import logsumexp
from numba import njit
R=Path('/mnt/data/egfr_campaign'); S=R/'stage2'; O=R/'intermediate'
sys.path.insert(0,str(R/'code'))
from legacy_geometry import atom_frame,unit
AA='ACDEFGHIKLMNPQRSTVWY'; BB=['N','CA','C','O']
LIB=json.loads((O/'empirical_rotamers.json').read_text())
RT=0.0019872041*298.15
ELEMENT_RADII={'C':1.70,'N':1.55,'O':1.52,'S':1.80}
BONDS={
 'A':[], 'G':[], 'C':[('CB','SG')], 'S':[('CB','OG')],
 'T':[('CB','OG1'),('CB','CG2')], 'V':[('CB','CG1'),('CB','CG2')],
 'I':[('CB','CG1'),('CB','CG2'),('CG1','CD1')],
 'L':[('CB','CG'),('CG','CD1'),('CG','CD2')],
 'M':[('CB','CG'),('CG','SD'),('SD','CE')],
 'D':[('CB','CG'),('CG','OD1'),('CG','OD2')],
 'N':[('CB','CG'),('CG','OD1'),('CG','ND2')],
 'E':[('CB','CG'),('CG','CD'),('CD','OE1'),('CD','OE2')],
 'Q':[('CB','CG'),('CG','CD'),('CD','OE1'),('CD','NE2')],
 'K':[('CB','CG'),('CG','CD'),('CD','CE'),('CE','NZ')],
 'R':[('CB','CG'),('CG','CD'),('CD','NE'),('NE','CZ'),('CZ','NH1'),('CZ','NH2')],
 'H':[('CB','CG'),('CG','ND1'),('CG','CD2'),('ND1','CE1'),('CD2','NE2'),('CE1','NE2')],
 'F':[('CB','CG'),('CG','CD1'),('CG','CD2'),('CD1','CE1'),('CD2','CE2'),('CE1','CZ'),('CE2','CZ')],
 'Y':[('CB','CG'),('CG','CD1'),('CG','CD2'),('CD1','CE1'),('CD2','CE2'),('CE1','CZ'),('CE2','CZ'),('CZ','OH')],
 'P':[('CB','CG'),('CG','CD'),('CD','N')],
 'W':[('CB','CG'),('CG','CD1'),('CG','CD2'),('CD1','NE1'),('NE1','CE2'),('CD2','CE2'),('CD2','CE3'),('CE2','CZ2'),('CE3','CZ3'),('CZ2','CH2'),('CZ3','CH2')]
}
CHARGE_NAMES={'D':['OD1','OD2'],'E':['OE1','OE2'],'K':['NZ'],'R':['NE','NH1','NH2'],'H':['ND1','NE2']}
CHARGES={'D':-1.,'E':-1.,'K':1.,'R':1.,'H':0.}
DONORS={'N','NE','NZ','NH1','NH2','ND2','NE2','NE1','OG','OG1','OH','ND1'}

def atom_data(aa,atoms,names=None):
    """xyz, element, radius, donor, acceptor, outward vector. No hydrogens."""
    if names is None:names=list(atoms)
    adjacency={n:[] for n in atoms}
    for a,b in [('N','CA'),('CA','C'),('C','O'),('CA','CB')]+BONDS[aa]:
        if a in atoms and b in atoms:adjacency[a].append(b);adjacency[b].append(a)
    ar=[]
    for n in names:
        p=np.array(atoms[n],float);e=n[0];don=n in DONORS
        if aa=='P' and n=='N':don=False
        acc=e=='O' or (aa=='H' and n in ['ND1','NE2']) or (aa=='C' and n=='SG')
        if aa=='Q' and n=='NE2':acc=False
        if aa=='R' and n=='NE':acc=False
        if aa=='H' and n in ['ND1','NE2']:don=True # neutral tautomer envelope; exact FF checked later
        nn=adjacency.get(n,[]);out=unit(p-np.mean([atoms[q] for q in nn],axis=0)) if nn else np.zeros(3)
        ar.append([*p,ELEMENT_RADII.get(e,1.7),float(e=='C'),float(don),float(acc),*out])
    return np.array(ar,dtype=float).reshape(-1,10)

def charge_center(aa,atoms):
    names=CHARGE_NAMES.get(aa,[])
    return np.mean([atoms[n] for n in names],axis=0) if names and all(n in atoms for n in names) else None

def conformer(aa,index,bb):
    names=LIB[aa]['names'];frame=atom_frame(*[np.array(bb[k]) for k in range(3)])
    xyz=np.array(LIB[aa]['xyz'][index],float).reshape(-1,3)@frame+bb[1]
    atoms={k:np.array(v) for k,v in zip(BB,bb[:4])}
    atoms.update({k:p for k,p in zip(names,xyz)})
    return atoms

@njit(cache=True)
def pair_score(a,b):
    """Returns steric penalty, contact attraction, directional hydrogen bonds."""
    rep=0.;vdw=0.;hb=0.
    for x in range(len(a)):
        for y in range(len(b)):
            dx=b[y,0]-a[x,0];dy=b[y,1]-a[x,1];dz=b[y,2]-a[x,2]
            d2=dx*dx+dy*dy+dz*dz
            if d2>42.25:continue
            d=math.sqrt(max(d2,0.04));r=a[x,3]+b[y,3]
            hard=.84*r
            if d<hard:rep+=35.*(hard-d)**2
            if d<1.6:rep+=500.*(1.6-d)**2+100.
            vdw-= (.10 if a[x,4]*b[y,4] else .035)*math.exp(-((d-r)/.75)**2)
            if d>2.1 and d<3.6:
                for flip in range(2):
                    if flip==0:don=a[x,5];acc=b[y,6];ux=a[x,7];uy=a[x,8];uz=a[x,9];vx=b[y,7];vy=b[y,8];vz=b[y,9];sgn=1.
                    else:don=b[y,5];acc=a[x,6];ux=b[y,7];uy=b[y,8];uz=b[y,9];vx=a[x,7];vy=a[x,8];vz=a[x,9];sgn=-1.
                    if don and acc:
                        c1=(ux*dx+uy*dy+uz*dz)*sgn/d;c2=-(vx*dx+vy*dy+vz*dz)*sgn/d
                        orient=max(0.,min(1.,(c1+.1)/.8))*max(0.,min(1.,(c2+.2)/.7))
                        hb-=.9*math.exp(-((d-2.85)/.42)**2)*orient*.5
    return rep,vdw,hb

@njit(cache=True)
def screened_charge(q1,q2,d,epsilon=20.):
    return 332.06371*q1*q2*math.exp(-.10*max(d-3.,0.))/(epsilon*max(d,2.5))

def all_atoms(rows):
    return np.concatenate([atom_data(r['aa'],r['atoms']) for r in rows])

def complete_framework(bb,mask,aa,orig,frames):
    """Repair incomplete framework side chains, retaining measured complete ones."""
    by={i:{} for i in range(len(bb))}
    for i,n,p in frames:by[i][n]=np.array(p)
    fixed={}
    for i in np.flatnonzero(~mask):
        if set(LIB[aa[i]]['names'])<=set(by[i]):fixed[i]=by[i]
    context=[atom_data(aa[i],a) for i,a in fixed.items()]
    for i in np.flatnonzero(~mask):
        if i in fixed:continue
        candidates=[]
        env=np.concatenate([atom_data(aa[j],a) for j,a in fixed.items() if abs(i-j)>1]+[atom_data('A',{n:p for n,p in zip(BB,bb[j])}) for j in np.flatnonzero(mask) if abs(i-j)>1])
        for k in range(len(LIB[aa[i]]['xyz'])):
            a=conformer(aa[i],k,bb[i]);dat=atom_data(aa[i],a,LIB[aa[i]]['names'])
            rep,vdw,hb=pair_score(dat,env)
            observed=sum(np.sum((a[n]-p)**2) for n,p in by[i].items() if n in a and n not in BB)
            candidates.append((rep+vdw+hb+observed,a,k))
        score,a,k=min(candidates,key=lambda t:t[0]);fixed[i]=a
    return fixed

def prepare_pose(rec):
    import dock_mt as dock
    bb,mask,aa,orig,pp,frames,rama,closure=dock.assemble(rec['selection'])
    rot=np.array(rec['rotation']);tr=np.array(rec['translation']);bb=bb@rot+tr
    frames=[(i,n,np.array(p)@rot+tr) for i,n,p in frames]
    fixed=complete_framework(bb,mask,aa,orig,frames)
    return bb,mask,aa,orig,pp,fixed

def proton_free_energy(shift,pH,pka=6.3):
    """Independent-site linkage; zero shift gives exactly zero binding effect."""
    x=math.log(10.)*(pka-pH)
    return -RT*(np.logaddexp(0.,x-np.array(shift)/RT)-np.logaddexp(0.,x))

@njit(cache=True)
def proton_energy_numba(shift,pH,pka=6.3):
    x=math.log(10.)*(pka-pH)
    a=x-shift/RT
    log1=(max(a,0.)+math.log1p(math.exp(-abs(a))))
    log0=(max(x,0.)+math.log1p(math.exp(-abs(x))))
    return -RT*(log1-log0)

@njit(cache=True)
def objective(sel,one,pair,ib,hs,tf,aaidx,pka,mode):
    """Joint human acid/mouse acid/neutral-human specificity surrogate.
    Independent-site calculation is a search approximation. Coupled enumeration
    is used for the later evaluation of actual sequences.
    """
    fold=0.;raw=np.zeros(2);th=np.zeros((2,tf.shape[3]));bindh=np.zeros(2)
    his=0;acid=0;aromatic=0;gly=0
    for i in range(len(sel)):
        oi=sel[i];fold+=one[i,oi]
        for j in range(i):fold+=pair[i,oi,j,sel[j]]
        ai=aaidx[i,oi]
        his+=int(ai==6);acid+=int(ai==2 or ai==3);aromatic+=int(ai==4 or ai==18 or ai==19);gly+=int(ai==5)
        for s in range(2):
            raw[s]+=ib[i,oi,s]
            if ai==6:bindh[s]+=proton_energy_numba(hs[i,oi,s],6.5,pka)
            for h in range(tf.shape[3]):th[s,h]+=tf[i,oi,s,h]
    low=raw.copy()+bindh
    high=raw[0]
    for i in range(len(sel)):
        if aaidx[i,sel[i]]==6:high+=proton_energy_numba(hs[i,sel[i],0],7.4,pka)
    for s in range(2):
        for h in range(th.shape[1]):low[s]+=proton_energy_numba(th[s,h],6.5,pka)
    for h in range(th.shape[1]):high+=proton_energy_numba(th[0,h],7.4,pka)
    switch=high-low[0]
    # Keep concentrations/affinity out: this is only a dimensioned heuristic.
    worst=max(low[0],low[1]); mismatch=abs(low[0]-low[1])
    composition=.6*max(his-6,0)**2+.28*max(acid-8,0)**2+.3*max(aromatic-10,0)**2+.12*max(gly-12,0)**2
    if mode==0:ans=fold+.80*worst+.30*mismatch-2.5*switch
    elif mode==1:ans=fold+.60*worst+.45*mismatch-4.*switch
    else:ans=fold+1.10*worst+.50*mismatch-1.8*switch
    return ans+composition,fold,low[0],low[1],high,switch,his

@njit(cache=True)
def anneal(nopt,one,pair,ib,hs,tf,aaidx,seed,steps=4000,mode=0,pka=6.3):
    # Numba's seeded np.random implementation is the legacy MT19937 generator.
    np.random.seed(seed);sel=np.zeros(len(nopt),dtype=np.int64)
    for i in range(len(sel)):sel[i]=np.random.randint(nopt[i])
    value=objective(sel,one,pair,ib,hs,tf,aaidx,pka,mode)[0];best=value;bs=sel.copy()
    for step in range(steps):
        i=np.random.randint(len(sel));old=sel[i];new=np.random.randint(nopt[i]);sel[i]=new
        nv=objective(sel,one,pair,ib,hs,tf,aaidx,pka,mode)[0]
        temp=2.8*(.08/2.8)**(step/max(steps-1,1))
        if nv<value or np.random.random()<math.exp(min(0.,(value-nv)/temp)):
            value=nv
            if value<best:best=value;bs=sel.copy()
        else:sel[i]=old
    # Deterministic one-site descent to eliminate obvious packing conflicts.
    sel=bs.copy()
    for sweep in range(3):
        changes=0
        for i in range(len(sel)):
            old=sel[i];v=objective(sel,one,pair,ib,hs,tf,aaidx,pka,mode)[0];winner=old
            for j in range(nopt[i]):
                sel[i]=j;nv=objective(sel,one,pair,ib,hs,tf,aaidx,pka,mode)[0]
                if nv<v:v=nv;winner=j
            sel[i]=winner;changes+=int(old!=winner)
        if not changes:break
    return sel,objective(sel,one,pair,ib,hs,tf,aaidx,pka,mode)
