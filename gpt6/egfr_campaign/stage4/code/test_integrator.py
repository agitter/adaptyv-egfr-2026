from pathlib import Path
import sys,numpy as np,json,math
S=Path('/mnt/data/egfr_campaign/stage4');sys.path.insert(0,str(S/'code'));from dynamics import mm,u,velocity_verlet
out=[]
sys1=mm.System();sys1.addParticle(12.);ff=mm.CustomExternalForce('500*(x*x+y*y+z*z)');ff.addParticle(0,[]);sys1.addForce(ff);it=velocity_verlet(.002);ctx=mm.Context(sys1,it,mm.Platform.getPlatformByName('Reference'));ctx.setPositions([[.1,0,0]]);ctx.setVelocities([[0,0,0]]);es=[]
for i in range(500):
 it.step(1);st=ctx.getState(getEnergy=True);es.append((st.getPotentialEnergy()+st.getKineticEnergy()).value_in_unit(u.kilojoule_per_mole))
drift=(max(es)-min(es))/np.mean(es);out.append({'test':'harmonic oscillator relative energy range','value':float(drift),'pass':bool(drift<.0001)});del ctx,it
sys2=mm.System();sys2.addParticle(12.);sys2.addParticle(1.);sys2.addConstraint(0,1,.1);it=velocity_verlet(.002);ctx=mm.Context(sys2,it,mm.Platform.getPlatformByName('Reference'));ctx.setPositions([[0,0,0],[.1,0,0]]);ctx.setVelocities([[0,.2,0],[0,-.2,0]]);ctx.applyVelocityConstraints(1e-8);it.step(1000);st=ctx.getState(getPositions=True,getVelocities=True);x=np.asarray(st.getPositions(asNumpy=True).value_in_unit(u.nanometer));v=np.asarray(st.getVelocities(asNumpy=True).value_in_unit(u.nanometer/u.picosecond));err=abs(np.linalg.norm(x[1]-x[0])-.1);dot=abs(float(np.dot(x[1]-x[0],v[1]-v[0])));out.extend([{'test':'RATTLE bond constraint','value':float(err),'pass':bool(err<1e-6)},{'test':'RATTLE tangent velocity constraint','value':dot,'pass':bool(dot<1e-6)}]);del ctx,it
(S/'reference/integrator_tests.json').write_text(json.dumps(out,indent=2));print(out);assert all(r['pass'] for r in out)
