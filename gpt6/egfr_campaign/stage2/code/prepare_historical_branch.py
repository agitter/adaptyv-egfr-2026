from pathlib import Path
import hashlib,json
R=Path('/mnt/data/egfr_campaign');S=R/'stage2'
s=(R/'code/03_loops.py').read_text()
s=s.replace("out=O/'loops'","out=R/'stage2/intermediate/loops'")
s=s.replace('rng=np.random.default_rng(20260929)','rng=np.random.Generator(np.random.MT19937(20031998))')
s=s.replace("R/'reference/loop_generation_stats.json'","R/'stage2/reference/loop_generation_stats.json'")
s=s.replace("(O/'framework_only.json').write_text(json.dumps(framework))","assert json.loads((O/'framework_only.json').read_text()) == framework")
s=s.replace("(R/'reference/framework_definition.json').write_text", "(R/'stage2/reference/framework_definition.json').write_text")
s=s.replace("np.save(O/'nonbinder_torsion_fragments.npy',frags)","np.save(R/'stage2/intermediate/nonbinder_torsion_fragments.npy',frags)")
(S/'code/generate_loops_mt.py').write_text(s)
s=(R/'code/05_dock.py').read_text()
s=s.replace("OUT=O/'docking'","OUT=R/'stage2/intermediate/docking'")
s=s.replace("(O/'loops').glob", "(R/'stage2/intermediate/loops').glob")
s=s.replace("R/'reference/docking_stats.json'", "R/'stage2/reference/docking_stats.json'")
# The previous clash filter omitted loop-to-nonadjacent framework contacts.
s=s.replace("where=mask[i]&mask[j]","where=mask[i]|mask[j]")
s=s.replace("others=mut[a+1:];others=others[others-i>1]","others=np.arange(len(bb));others=others[(np.abs(others-i)>1)&((~mask[others])|(others>i))]")
s=s.replace("if dd.min()<2.15", "if dd.min()<2.15")
(S/'code/dock_mt.py').write_text(s)
(S/'reference/historical_branch.json').write_text(json.dumps({'seed_loop':20031998,'seed_docking':19762003,'rng':'MT19937 (1998)','ancestry':'New loops and docking; earlier PCG64 pilot excluded','added_filter':'all loop versus nonadjacent framework backbone pairs, not only loop-loop','sources':{'loop_parent_sha256':hashlib.sha256((R/'code/03_loops.py').read_bytes()).hexdigest(),'docking_parent_sha256':hashlib.sha256((R/'code/05_dock.py').read_bytes()).hexdigest()}},indent=2))
