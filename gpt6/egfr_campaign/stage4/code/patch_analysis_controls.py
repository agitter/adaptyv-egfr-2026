"""One-time auditable analysis correction; does not modify trajectory data."""
from pathlib import Path
import shutil,json
S=Path('/mnt/data/egfr_campaign/stage4')
p=S/'code/trajectory_contacts.py'
shutil.copy2(p,S/'reference/trajectory_contacts_v2_before_matched_control.py')
t=p.read_text()
a="root=next(r for r in [S,R/'stage3'] if (r/'intermediate/designs'/(cid+'.json')).exists());design=json.loads((root/'intermediate/designs'/(cid+'.json')).read_text());"
b="actual=cid.removeprefix('MATCHED_');root=S/'intermediate/matched_controls' if cid.startswith('MATCHED_') else next(r for r in [S,R/'stage3'] if (r/'intermediate/designs'/(actual+'.json')).exists());design=json.loads((root/'intermediate/designs'/(actual+'.json')).read_text());"
assert a in t;t=t.replace(a,b);p.write_text(t)
p=S/'code/analyze_trajectories.py';shutil.copy2(p,S/'reference/analyze_trajectories_v1_reference_mismatch.py');t=p.read_text()
t=t.replace("qualification='Native control reference is the initial recorded, minimized frame, not pre-repair crystal coordinates.';geom=None", "qualification='Reference for all internal-deformation comparisons is the initial recorded, minimized unbound frame.';native_sequence=''.join(__import__('Bio.SeqUtils',fromlist=['seq1']).seq1(r.name) for r in __import__('trajectory_contacts').app.PDBFile(str(final)).topology.residues());geom=audit(final,{'sequence':native_sequence,'candidate_id':cid},chain='A')")
t=t.replace("assert len(intervals)==3;geom=audit(final,design,chain='A');qualification='Design reference is its separately refined complex geometry, before the unbound preparation minimization.'", "assert len(intervals)==3;geom=audit(final,design,chain='A');reference=coords[0];qualification='Reference for all internal-deformation comparisons is the initial recorded, minimized unbound frame.'")
t=t.replace("'end_point_independent_geometry':geom,", "'end_point_independent_geometry':geom,'end_point_geometry_qualification':'The omega/disulfide thresholds were designed for minimized structures. Thermal excursions at one instantaneous endpoint are diagnostic only, not automatic exclusion. Chirality inversion or severe bond/clash defects remain separate warnings.',")
p.write_text(t)
(S/'reference/trajectory_analysis_reference_correction.json').write_text(json.dumps({'purpose':'Use the same initial recorded unbound minimized reference for native and designed loop internal deformation; add matched-parent endpoint support. Qualify instantaneous thermal geometry warnings separately from minimized-structure gates.','raw_trajectory_files_modified':False,'superseded_code_preserved':True},indent=2))
