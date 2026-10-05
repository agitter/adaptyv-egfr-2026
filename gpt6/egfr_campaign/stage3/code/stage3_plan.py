"""Record the falsification-driven stage-3 plan and immutable input ancestry."""
from pathlib import Path
import json, hashlib, time
R=Path('/mnt/data/egfr_campaign');S=R/'stage3'
plan={
 'date':'2026-10-01','stage':'3','status':'in progress','phone_a_friend_used':2,
 'goals':['Improve actual sequence/structure candidates; no invented affinity or eligibility claims.',
 'Receptor-histidine recognition, with both human and mouse favorable states but human neutral-pH negative design.',
 'Repair initial search omissions before producing another large library.',
 'Evaluate several redesigned interfaces with matched-component Amber99SB/OBC, protonation and neutral-tautomer sensitivity.'],
 'corrections':['Explicit intrabinder fixed-charge electrostatics in one- and two-body packing terms.',
 'Target-attached glycans and complete fixed framework atoms checked during rotamer selection.',
 'Asymmetric mouse objective: stronger mouse recognition is not penalized.',
 'Explicit interface desolvation and restrained charge density; no reward for histidine abundance.',
 'Geometry/conformation, matched protonation energies and unbound normalization remain separate gates.'],
 'scientific_methods':[{'name':'Coulomb/Debye electrostatics and dielectric sensitivity','not_after':1923},
 {'name':'Empirical rotamers, rigid docking, Metropolis/annealing and fixed-backbone design','not_after':2000},
 {'name':'Cyclic coordinate descent loop closure','not_after':2003},
 {'name':'Shrake-Rupley surface accessibility','not_after':1973},
 {'name':'MT19937','not_after':1998},
 {'name':'Amber99SB / OBC generalized Born','not_after':2006},
 {'name':'Binding polynomial / Wyman proton linkage','not_after':1964},
 {'name':'Discrete protonation-state GB calculations','not_after':2004},
 {'name':'BLAST and Levenshtein novelty screens','not_after':1990}],
 'qualification':'Modern implementations and supplied modern data are permitted. No learned folding or protein generator is used. Earlier excluded PCG64 pilot remains disclosed.',
 'sources':[
 {'url':'https://proteinbase.com/competitions/anthropic-adaptyv-2026/challenges/egfr','role':'live fetch unavailable this stage; supplied HTML and FAQ used'},
 {'url':'https://www.adaptyvbio.com/blog/novelty','role':'current novelty rules, retrieved 2026-10-01'},
 {'url':'https://doi.org/10.1002/jcc.20139','role':'Mongan/Case/McCammon discrete protonation with generalized Born, 2004'},
 {'url':'https://doi.org/10.1016/S0065-3233(08)60190-4','role':'Wyman linkage, 1964'},
 {'url':'https://doi.org/10.1016/S1097-2765(01)00230-1','role':'FcRn control; 2001'},
 {'url':'https://doi.org/10.1038/nsb877','role':'Havranek and Harbury 2003 explicit positive/negative design, verified'}],
 'reference_files':[]}
for name in ['stage2/STATE.md','stage2/code/design_search.py','stage2/code/classical_design.py','stage2/code/mm_refine.py','stage2/code/mm_binding_probe.py','stage2/intermediate/docking/poses.json','inputs/provided/challenge-faq.md','intermediate/challenge_saved.txt']:
 p=R/name;plan['reference_files'].append({'path':name,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
(S/'reference/plan.json').write_text(json.dumps(plan,indent=2))
