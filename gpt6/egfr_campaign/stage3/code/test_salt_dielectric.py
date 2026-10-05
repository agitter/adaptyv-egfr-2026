"""Verify the classical salt/dielectric extension against independent OpenMM."""
from fast_obc import *
from openmm.app.internal.customgbforces import GBSAOBC2Force

def engine(x,q,sigma,eps,radius,scale,dielectric,kappa):
    system=mm.System();nb=mm.NonbondedForce();gb=GBSAOBC2Force(solventDielectric=78.3,soluteDielectric=dielectric,SA='ACE',kappa=kappa)
    for i in range(len(x)):
        system.addParticle(12.);nb.addParticle(q[i]/math.sqrt(dielectric),sigma[i],eps[i]);gb.addParticle([q[i],radius[i],scale[i]])
    gb.finalize();system.addForce(nb);system.addForce(gb);integ=mm.VerletIntegrator(.001);ctx=mm.Context(system,integ,mm.Platform.getPlatformByName('Reference'));ctx.setPositions(x)
    val=float(ctx.getState(getEnergy=True).getPotentialEnergy().value_in_unit(u.kilojoule_per_mole));del ctx,integ;return val
rng=np.random.Generator(np.random.MT19937(20041998));x=np.array([[i*.37,j*.43,k*.45] for i in range(3) for j in range(2) for k in range(2)]);x+=rng.uniform(-.03,.03,x.shape)
q=rng.uniform(-.8,.8,len(x));sigma=np.full(len(x),.3);eps=np.full(len(x),.2);radius=np.full(len(x),.16);scale=np.full(len(x),.8);part=np.array([0]*6+[1]*6,np.int64)
lj,sa,A,B,*_=interaction_coefficients(x,q,sigma,eps,radius,scale,part,np.array([0.,1.25]));tests=[]
for de in [1.,2.,4.]:
    for ki,kap in enumerate([0.,1.25]):
        native=engine(x,q,sigma,eps,radius,scale,de,kap)
        for p in [0,1]:
            mask=part==p;native-=engine(x[mask],q[mask],sigma[mask],eps[mask],radius[mask],scale[mask],de,kap)
        direct=lj+sa+A/de+B[ki]/78.3;error=abs(native-direct)/4.184
        tests.append({'dielectric':de,'kappa':kap,'error_kcal':float(error),'pass':bool(error<1e-4)})
record={'tests':tests,'all_pass':all(t['pass'] for t in tests)};(S3/'reference/salt_dielectric_tests.json').write_text(json.dumps(record,indent=2));print(record);assert record['all_pass']
