"""Classical internal coordinates, Kabsch fitting, and CCD loop closure.
Scientific algorithms: Kabsch 1976; Canutescu & Dunbrack 2003.
Numba is used only as a modern compiler, not as a scientific model.
"""
import numpy as np, math
from numba import njit

@njit(cache=True)
def unit(v):
    return v / max(np.sqrt(np.dot(v,v)),1e-12)

@njit(cache=True)
def dihedral(a,b,c,d):
    b0=a-b; b1=unit(c-b); b2=d-c
    v=b0-np.dot(b0,b1)*b1;w=b2-np.dot(b2,b1)*b1
    return math.atan2(np.dot(np.cross(b1,v),w),np.dot(v,w))

@njit(cache=True)
def place(a,b,c,length,angle,torsion):
    bc=unit(c-b); normal=unit(np.cross(b-a,bc));m=np.cross(normal,bc)
    return c+length*(-bc*math.cos(angle)+m*math.sin(angle)*math.cos(torsion)+normal*math.sin(angle)*math.sin(torsion))

@njit(cache=True)
def make_chain(anchor,phis,psis,psi_anchor,target):
    n=len(phis)-1; x=np.zeros(((n+2)*3,3));x[:3]=anchor
    for j in range(n+1):
        k=3*(j+1); psi=psi_anchor if j==0 else psis[j-1]
        x[k]=place(x[k-3],x[k-2],x[k-1],1.329,math.radians(116.2),psi)
        l1=1.458;l2=1.525;ang=math.radians(111.2)
        if j==n:
            l1=np.linalg.norm(target[1]-target[0]);l2=np.linalg.norm(target[2]-target[1]);ang=math.acos(np.dot(unit(target[0]-target[1]),unit(target[2]-target[1])))
        x[k+1]=place(x[k-2],x[k-1],x[k],l1,math.radians(121.7),math.pi)
        x[k+2]=place(x[k-1],x[k],x[k+1],l2,ang,phis[j])
    return x

@njit(cache=True)
def ccd(x,target,max_cycles=450,tolerance=0.075):
    n=len(x)//3-2
    for cycle in range(max_cycles):
        for res in range(1,n+2):
            for off in range(2):
                k=res*3+off
                if k+2>=len(x):continue
                a=x[k].copy();u=unit(x[k+1]-a);num=0.;den=0.
                for e in range(3):
                    v=x[len(x)-3+e]-a; w=target[e]-a
                    vp=v-np.dot(v,u)*u;wp=w-np.dot(w,u)*u
                    num+=np.dot(u,np.cross(vp,wp));den+=np.dot(vp,wp)
                angle=math.atan2(num,den)
                angle=max(-0.25,min(0.25,angle))
                ca=math.cos(angle);sa=math.sin(angle)
                for j in range(k+2,len(x)):
                    v=x[j]-a
                    x[j]=a+v*ca+np.cross(u,v)*sa+u*np.dot(u,v)*(1-ca)
        err=np.sqrt(np.sum((x[-3:]-target)**2)/3)
        if err<tolerance:return x,err,cycle+1
    return x,err,max_cycles

@njit(cache=True)
def backbone_atoms(x):
    n=len(x)//3-2; out=np.zeros((n,5,3)); pp=np.zeros((n,2))
    for i in range(n):
        j=3*(i+1);N=x[j];CA=x[j+1];C=x[j+2];Nnext=x[j+3]
        phi=dihedral(x[j-1],N,CA,C);psi=dihedral(N,CA,C,Nnext)
        out[i,0]=N;out[i,1]=CA;out[i,2]=C
        out[i,3]=place(N,CA,C,1.231,math.radians(120.8),psi+math.pi)
        u=unit(N-CA);v=unit(C-CA);bis=unit(u+v);norm=unit(np.cross(u,v))
        out[i,4]=CA+1.53*(-0.588*bis+0.809*norm)
        pp[i,0]=phi;pp[i,1]=psi
    return out,pp

@njit(cache=True)
def rama_quality(pp):
    centers=np.array([[-65.,-40.,22.,25.],[-120.,135.,35.,30.],[-70.,145.,22.,30.],[65.,40.,25.,35.],[80.,-25.,25.,30.],[-80.,-170.,25.,25.]])
    out=np.zeros(len(pp));positive=0
    for i in range(len(pp)):
        a=pp[i,0]*180/math.pi;b=pp[i,1]*180/math.pi
        if a>0:positive+=1
        best=1000.
        for k in range(len(centers)):
            da=(a-centers[k,0]+180)%360-180;db=(b-centers[k,1]+180)%360-180
            d=(da/centers[k,2])**2+(db/centers[k,3])**2
            best=min(best,d)
        out[i]=best
    return out,positive

@njit(cache=True)
def loop_clashes(bb,framework,fids,left,right):
    bad=0;mind=99.
    for i in range(len(bb)):
        for a in range(4):
            for k in range(len(framework)):
                if i<2 and fids[k]==left:continue
                if i>=len(bb)-2 and fids[k]==right:continue
                d=np.linalg.norm(bb[i,a]-framework[k]);mind=min(mind,d)
                if d<2.15:bad+=1
            for j in range(i+2,len(bb)):
                for b in range(4):
                    d=np.linalg.norm(bb[i,a]-bb[j,b]);mind=min(mind,d)
                    if d<2.15:bad+=1
    return bad,mind

def kabsch(a,b):
    u,s,v=np.linalg.svd((a-a.mean(0)).T@(b-b.mean(0)))
    m=np.eye(3);m[2,2]=np.linalg.det(u@v);rot=u@m@v
    return rot,b.mean(0)-a.mean(0)@rot

def atom_frame(n,ca,c):
    x=unit(c-ca);z=unit(np.cross(n-ca,c-ca));y=np.cross(z,x)
    return np.array([x,y,z])
