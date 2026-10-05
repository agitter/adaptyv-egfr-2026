"""Write the stage's source-use and historical-method audit records."""
from pathlib import Path
import json,datetime
S=Path('/mnt/data/egfr_campaign/stage4')
sources=[
 {'id':'competition','year':2026,'url':'https://proteinbase.com/competitions/anthropic-adaptyv-2026/challenges/egfr','role':'Current rules, not a scientific algorithm','findings':['Human binding at pH6.5 and no detectable binding at pH7.4; mouse cross-reactivity; human affinity.','Track3 max20; single-chain10-250aa; nanobody common framework allowed.','No existing binder seed; CDR novelty checked; submission CSV name,sequence,molecule_class.','Deadline October4 2026 23:59AoE. No numeric binding-detection threshold specified on page.']},
 {'id':'novelty','year':2026,'url':'https://www.adaptyvbio.com/blog/novelty','role':'Organizer classification/screening; not reproduced exactly','findings':['Novel antibody CDR3 with familiar framework is treated separately from generic-protein novelty.','Exact official annotation/database snapshots not locally replicated.']},
 {'id':'entropy_stabilization','year':1987,'doi':'10.1073/pnas.84.19.6663','url':'https://pmc.ncbi.nlm.nih.gov/articles/PMC299143/','role':'Rationale for geometry-compatible Gly-to-Ala and X-to-Pro substitutions','qualification':'Published lysozyme results do not quantify effects in these designed nanobodies.'},
 {'id':'nanobody_disulfides','year':2008,'doi':'10.1016/j.jmb.2008.01.022','url':'https://europepmc.org/article/MED/18262543','role':'Rationale and caution for extra cystines','qualification':'Only one of two tested crosslink sites consistently stabilized the published nanobodies. Our57-60site is a different, unvalidated hypothesis.'},
 {'id':'canonical_loops','year':1987,'doi':'10.1016/0022-2836(87)90412-8','url':'https://www.sciencedirect.com/science/article/abs/pii/0022283687904128','role':'Considered next strategy; NOT implemented in stage4','qualification':'Limited canonical classes in five hypervariable loops; sequence/packing/torsion constraints matter, not loop length alone.'},
 {'id':'wyman_linkage','year':1964,'doi':'10.1016/S0065-3233(08)60190-4','url':'https://www.sciencedirect.com/science/chapter/bookseries/pii/S0065323308601904','role':'Binding-polynomial/proton-linkage basis'},
 {'id':'linked_protonation','year':1996,'doi':'10.1016/S0006-3495(96)79403-1','url':'https://pmc.ncbi.nlm.nih.gov/articles/PMC1233671/','role':'Classical proton-linked binding and intrinsic versus observed energetics','qualification':'We did not run calorimetry; no experimental pKa or affinity was inferred.'},
 {'id':'andersen','year':1980,'doi':'10.1063/1.439486','role':'Thermal collisions; implemented as constraint-group velocity collisions'},
 {'id':'rattle','year':1983,'doi':'10.1016/0021-9991(83)90014-1','role':'Synchronized constrained dynamics'},
 {'id':'leapfrog_docs','url':'https://docs.openmm.org/latest/api-python/generated/openmm.openmm.VerletIntegrator.html','role':'Modern implementation documentation; half-step velocity timing explains excluded pilot bias'},
 {'id':'ccd','year':2003,'doi':'10.1110/ps.0242703','role':'Inherited de novo loop construction'},
 {'id':'multistate','year':2003,'doi':'10.1038/nsb877','role':'Inherited positive/negative multistate sequence design'},
 {'id':'amber99sb','year':2006,'doi':'10.1002/prot.21123','role':'Forcefield'},
 {'id':'obc','year':2004,'doi':'10.1002/prot.20033','role':'Implicit solvent; Onufriev,Bashford,Case2004'},
 {'id':'framework','year':2008,'doi':'10.1074/jbc.M806889200','url':'https://www.rcsb.org/structure/3EAK','role':'Framework-only scaffold; native loops only in validation control'},
]
(S/'reference/sources.json').write_text(json.dumps({'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'source_use':'Primary-source abstract/HTML reading; no newly downloaded PDF or learned protein service. DOIs from inherited audit are preserved as scientific method identifiers.','sources':sources},indent=2))
record={'priority_order':['human acidic-selective binding/no detectable neutral binding','mouse cross-reactivity','human affinity'],'format':{'min_length':10,'max_length':250,'track3_max_upload':20,'requested_ranked_pool':100,'actual_submission_format':'CSV(name,sequence,molecule_class)','requested_deliverable':'FASTA plus compatible CSV; pending'},'all_stage4_sequences_length':125,'common_framework_allowed':True,'no_existing_binder_seed':True,'current_designs_derive_from_own_de_novo_parents':True,'local_novelty_not_official_acceptance':True,'biological_criteria_experimentally_established':False,'neutral_no_detection_established':False,'phone_a_friend_used':2,'phone_a_friend_remaining':1}
(S/'reference/criteria_status.json').write_text(json.dumps(record,indent=2))
print('Wrote',len(sources),'source records and criteria status')
