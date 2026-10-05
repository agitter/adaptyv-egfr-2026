from pathlib import Path
import json
S=Path('/mnt/data/egfr_campaign/stage4');p=S/'code/dynamics.py';text=p.read_text();(S/'reference/dynamics_v1_leapfrog_temperature_bias.py').write_text(text)
func='''def velocity_verlet(dt_ps):
    """Synchronized velocity Verlet with RATTLE-style constraint corrections.
    Scientific equations predate2011; CustomIntegrator is only the implementation.
    """
    it=mm.CustomIntegrator(dt_ps*u.picosecond)
    it.addPerDofVariable('x_free',0.)
    it.addUpdateContextState()
    it.addComputePerDof('v','v+0.5*dt*f/m')
    it.addComputePerDof('x','x+dt*v')
    it.addComputePerDof('x_free','x')
    it.addConstrainPositions()
    it.addComputePerDof('v','v+0.5*dt*f/m+(x-x_free)/dt')
    it.addConstrainVelocities()
    it.setConstraintTolerance(1e-6)
    return it

'''
text=text.replace('def run(cid,seed=',func+'def run(cid,seed=')
text=text.replace("tag=f'{cid}_s{seed}'","tag=f'{cid}_v2_s{seed}'")
text=text.replace("intervals=[[25,36],[50,65],[97,116]]","intervals=[[25,38],[52,67],[99,117]]")
text=text.replace("'native index regions approximating design boundaries; loop lengths differ'","'native source residue regions26-38,53-67,100-117; loop lengths differ from designs'")
text=text.replace("integrator=mm.VerletIntegrator(dt*u.picosecond)","integrator=velocity_verlet(dt)")
text=text.replace("'protocol':'Andersen-style constraint-group collisions every0.05ps,5/ps; Verlet;300K;NoCutoff'","'protocol':'v2:Andersen-style constraint-group collisions every0.05ps,5/ps;synchronized velocity-Verlet/RATTLE;300K;NoCutoff'")
p.write_text(text)
(S/'reference/dynamics_protocol_correction.json').write_text(json.dumps({'old_protocol':'leapfrog stored half-step velocities were Maxwell refreshed without time synchronization','observed_pilot_production_temperatures_K':[309.3952,311.3236,311.9632],'fix':'Use synchronized velocity-Verlet plus RATTLE-style projection with the same collision rule. Correct native CDR masks to the historical framework deletion regions.','historical_methods':['velocity Verlet pre1983','RATTLE1983','Andersen1980','MT199371998'],'pilot_outputs':'retained, excluded from corrected comparison'},indent=2))
