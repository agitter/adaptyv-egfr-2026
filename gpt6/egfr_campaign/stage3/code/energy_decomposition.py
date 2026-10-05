"""Exact additive accounting of the implemented matched-coordinate OBC proxy.
Half of each pair term is assigned to each residue; this accounting convention
is not a causal per-residue binding free energy. It guides mutation hypotheses.
"""
from fast_obc import *

@njit(cache=True)
def decompose(xyz,q,sigma,eps,radius,scale,partner,resid,kappas):
    bf,bo=born_pair(xyz,radius,scale,partner);n=len(q);out=np.zeros((resid.max()+1,3+len(kappas)))
    for i in range(n):
        ri=resid[i];out[ri,1]+=4*math.pi*2.25936*(radius[i]+.14)**2*((radius[i]/bf[i])**6-(radius[i]/bo[i])**6)
        out[ri,2]-=.5*K*q[i]*q[i]*(1/bf[i]-1/bo[i])
        for t in range(len(kappas)):out[ri,3+t]+=.5*K*q[i]*q[i]*(math.exp(-kappas[t]*bf[i])/bf[i]-math.exp(-kappas[t]*bo[i])/bo[i])
        for j in range(i):
            rj=resid[j];delta=xyz[i]-xyz[j];r2=np.dot(delta,delta);r=math.sqrt(r2);bb=bf[i]*bf[j];f=math.sqrt(r2+bb*math.exp(-r2/(4*bb)));c=K*q[i]*q[j];val=np.zeros(3+len(kappas))
            if partner[i]!=partner[j]:
                sr=.5*(sigma[i]+sigma[j])/r;sr6=sr**6;val[0]=4*math.sqrt(eps[i]*eps[j])*(sr6*sr6-sr6);val[2]=c/r-c/f
                for t in range(len(kappas)):val[3+t]=c*math.exp(-kappas[t]*f)/f
            else:
                b0=bo[i]*bo[j];f0=math.sqrt(r2+b0*math.exp(-r2/(4*b0)));val[2]=-c*(1/f-1/f0)
                for t in range(len(kappas)):val[3+t]=c*(math.exp(-kappas[t]*f)/f-math.exp(-kappas[t]*f0)/f0)
            out[ri]+=.5*val;out[rj]+=.5*val
    return out

def main(cid):
    src=S3/'intermediate/refined'/(cid+'_human6ARU_allHIP.pdb');p=app.PDBFile(str(src));pars=parameters(p.topology,p.positions);res=np.array([a.residue.index for a in p.topology.atoms()],np.int64);kap=np.array([0.,1.25]);arr=decompose(*pars,res,kap);lj,sa,A,B,*_=interaction_coefficients(*pars,kap);expected=np.r_[lj,sa,A,B];error=float(np.max(np.abs(arr.sum(0)-expected)));assert error<1e-6
    rows=[]
    for r,a in zip(p.topology.residues(),arr):
        rows.append({'chain':r.chain.id,'position':int(r.id),'residue':r.name,'LJ_kcal':float(a[0]/4.184),'SA_kcal':float(a[1]/4.184),'epsilon1_salt_kcal':float((a[0]+a[1]+a[2]+a[4]/78.3)/4.184),'epsilon4_salt_kcal':float((a[0]+a[1]+a[2]/4.+a[4]/78.3)/4.184)})
    rec={'source':str(src),'sum_check_max_error_kj':error,'interpretation':'Pairwise accounting only; mutation effects require recomputation. All selected histidines are HIP in this diagnostic.','residues':rows};out=S3/'intermediate/decomposition';out.mkdir(exist_ok=True);(out/(cid+'.json')).write_text(json.dumps(rec,indent=2))
    for r in sorted([r for r in rows if r['chain']=='B'],key=lambda r:r['epsilon1_salt_kcal'],reverse=True)[:16]:print(cid,r)
if __name__=='__main__':main(sys.argv[1])
