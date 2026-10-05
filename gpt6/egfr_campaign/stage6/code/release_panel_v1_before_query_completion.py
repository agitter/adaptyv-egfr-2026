"""Freeze, export, and independently validate the ranked experimental panel.
No measured affinity, pH selectivity, expression, or organizer eligibility is
inferred from this release. CSVs are machine interchange, not formatted sheets.
"""
from pathlib import Path
import csv, collections, hashlib, io, json, re, shutil, sys, textwrap, datetime
import numpy as np
from select_panel import make_pool, rank_pool
from provenance_audit import audit as ancestry_audit
R=Path('/mnt/data/egfr_campaign'); S=R/'stage6'; O=S/'output'
AA=set('ACDEFGHIKLMNPQRSTVWY')
def read(p): return json.loads(Path(p).read_text())
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def parse_fasta(path):
    rows=[]; name=None; seq=[]
    for line in Path(path).read_text().splitlines():
        if line.startswith('>'):
            if name is not None: rows.append((name,''.join(seq)))
            name=line[1:]; seq=[]
        elif line.strip():
            assert name is not None, 'Sequence precedes FASTA header'
            seq.append(line.strip())
    if name is not None: rows.append((name,''.join(seq)))
    return rows

def main():
    O.mkdir(exist_ok=True)
    batches={}
    for stem in ['human_factorial','mouse_factorial','cluster_acids','diversity_followup','matched_unbound','mouse_completion']:
        d=read(S/'reference'/f'{stem}_execution.json'); assert d['finished'] and d['completed']==d['jobs'],stem
        if stem!='matched_unbound': assert all(r['returncode']==0 for r in d['records']),(stem,d)
        else:
            recovery=read(S/'reference/unbound_chain_audit_recovery.json')
            assert recovery['recovered_checks']==10 and all(r['pass'] for r in recovery['results'])
            assert all(r['returncode']==1 for r in d['records'])
        batches[stem]={'jobs':d['jobs'],'finished':True,'original_returncodes':dict(collections.Counter(str(r['returncode']) for r in d['records']))}
    assert read(S/'reference/panel_blast_complete.json')['complete']
    pool,excluded=make_pool(); panel=rank_pool(pool,100); assert len(panel)==100
    # Every selected sequence needs an explicit glycan result, not an absent key.
    assert all(r['glycan_stress'] and r['glycan_stress']['conservative_panel_filter_pass'] for r in panel)
    # Ensure refresh/selection cannot accidentally resurrect old glycan failures.
    ungated,_=make_pool(apply_glycan_filter=False)
    gs={r['sequence_sha256']:r for r in read(S/'reference/glycan_pool_screen.json')['records']}
    assert all(r['sequence_sha256'] in gs for r in ungated)
    pool2,excluded2=make_pool(); panel2=rank_pool(pool2,100)
    assert [r['sequence_sha256'] for r in panel]==[r['sequence_sha256'] for r in panel2]
    search_rows=read(S/'reference/novelty_query_superset.json'); queries={r['sequence']:r['query_id'] for r in search_rows}
    missing=[r['candidate_key'] for r in panel if r['sequence'] not in queries]
    assert not missing,('Final panel outside BLAST superset',missing)
    dbs={n:read(S/'reference'/f'panel_blast_{n}_summary.json') for n in ['swissprot','pdb','antibodies_augmented']}
    for name,d in dbs.items():
        assert d['all_query_ids_verified'] and d['positive_control_exact_match'],name
    # Public codes disclose order, but not ancestry, epitope, or design strategy.
    stems=['AmideAlibi','FoldAndDagger','Covertase','PeptidePhantom','ResidueRuse',
           'IonCurtain','ChiralCipher','MolarMirage','EnzymeEnigma','SilentLigand',
           'BufferBandit','HelixHush','CovalentCover','KineticCloak','BondVoyage',
           'Stereosecret','IsoelectricAlias','CatalyticCaper','MotifMisdirect','PolypeptidePlot']
    covers=['Onyx','Velvet','Cipher','Cinder','Nightfall']
    aliases=[a+'_'+b for a in stems for b in covers]
    rng=np.random.RandomState(20061973); aliases=[aliases[int(i)] for i in rng.permutation(len(aliases))]
    for r,alias in zip(panel,aliases):
        r['public_name']=f'{r["rank"]:03d}_{alias}'
        assert re.fullmatch(r'[A-Za-z0-9_]+',r['public_name']) and len(r['public_name'])<=40
        assert r['length']==len(r['sequence']) and 10<=r['length']<=250 and set(r['sequence'])<=AA
        r['BLAST_query_id']=queries[r['sequence']]
        r['BLAST_high_coverage_hits']={n:[h for h in d['hits'][queries[r['sequence']]] if h['query_coverage_pct']>=90] for n,d in dbs.items()}
        r['ranking_scope']='Experimental priority, evidence-stratified. Within a stratum use conservative pH proxy with a declared CDR3-redundancy penalty; not a predicted competition outcome or measured binding quality.'
    assert len({r['sequence'] for r in panel})==len({r['public_name'] for r in panel})==100
    assert max(collections.Counter(r['CDRs'][2] for r in panel).values())<=6
    # No false claim that all three biological criteria have passed.
    assert all(all(v=='UNVERIFIED' for v in r['biological_criteria'].values()) for r in panel)
    anc=ancestry_audit(panel); assert anc['all_pass'] and anc['count']==100
    (S/'reference/final_ancestry_audit.json').write_text(json.dumps(anc,indent=2))
    (O/'final_panel_private.json').write_text(json.dumps(panel,indent=2))
    for filename,rows in [('egfr_ranked_100.fasta',panel),('egfr_track3_top20.fasta',panel[:20])]:
        (O/filename).write_text(''.join('>'+r['public_name']+'\n'+'\n'.join(textwrap.wrap(r['sequence'],80))+'\n' for r in rows))
    columns=['name','sequence','molecule_class']
    for filename,rows in [('egfr_track3_top20.csv',panel[:20]),('egfr_100_RESERVE_POOL_NOT_SINGLE_UPLOAD.csv',panel)]:
        with (O/filename).open('w',newline='') as f:
            writer=csv.DictWriter(f,fieldnames=columns,lineterminator='\n'); writer.writeheader()
            writer.writerows({'name':r['public_name'],'sequence':r['sequence'],'molecule_class':'nanobody'} for r in rows)
    fields=['rank','public_name','internal_key','sequence_sha256','length','evidence_tier','evidence_label','conservative_pH_proxy_kcal','human_passing_models','mouse_passing_models','local_CDR3_max_edit_identity','official_novelty_status','CDR1','CDR2','CDR3','glycan_unweighted_overlap_fraction','experimental_status']
    with (O/'private_codebook_and_evidence.tsv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n');w.writeheader()
        for r in panel:
            w.writerow({'rank':r['rank'],'public_name':r['public_name'],'internal_key':r['candidate_key'],'sequence_sha256':r['sequence_sha256'],'length':r['length'],'evidence_tier':r['tier'],'evidence_label':r['evidence_label'],'conservative_pH_proxy_kcal':f'{r["ranking_contrast_kcal_proxy"]:.6f}','human_passing_models':len(r['human_models']),'mouse_passing_models':len(r['mouse_models']),'local_CDR3_max_edit_identity':f'{r["local_novelty"]["max_edit_identity"]:.6f}','official_novelty_status':r['official_novelty_status'],'CDR1':r['CDRs'][0],'CDR2':r['CDRs'][1],'CDR3':r['CDRs'][2],'glycan_unweighted_overlap_fraction':f'{r["glycan_stress"]["worst_unweighted_overlap_fraction"]:.6f}','experimental_status':'UNVERIFIED: human acid binding, neutral nondetection, mouse binding, affinity, expression/folding'})
    # Independent format parsing and cross-file order/sequence checks.
    assert parse_fasta(O/'egfr_ranked_100.fasta')==[(r['public_name'],r['sequence']) for r in panel]
    assert parse_fasta(O/'egfr_track3_top20.fasta')==parse_fasta(O/'egfr_ranked_100.fasta')[:20]
    for name,count in [('egfr_track3_top20.csv',20),('egfr_100_RESERVE_POOL_NOT_SINGLE_UPLOAD.csv',100)]:
        with (O/name).open(newline='') as f:
            reader=csv.DictReader(f); assert reader.fieldnames==columns; rows=list(reader)
        assert len(rows)==count
        assert [(r['name'],r['sequence']) for r in rows]==parse_fasta(O/'egfr_ranked_100.fasta')[:count]
        assert all(r['molecule_class']=='nanobody' for r in rows)
    # Explicit novelty/sequence gates, while preserving the official status.
    assert all(r['local_novelty']['max_edit_identity']<.7 for r in panel)
    assert all(r['sequence'].count('C')==2 and not re.search(r'N[^P][ST]',r['sequence']) for r in panel)
    assert all(r['official_novelty_status']=='PENDING_ADAPTYV' for r in panel)
    frozen={p.name:{'bytes':p.stat().st_size,'sha256':sha(p)} for p in O.iterdir() if p.is_file() and p.suffix in ['.fasta','.csv','.tsv','.json'] and p.name!='release_validation.json'}
    report={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'count':100,'unique_sequences':100,'names_unique':100,'length_range':[min(r['length'] for r in panel),max(r['length'] for r in panel)],'top20_count':20,'molecule_class':'nanobody','evidence_tiers_all100':dict(collections.Counter(str(r['tier']) for r in panel)),'evidence_tiers_top20':dict(collections.Counter(str(r['tier']) for r in panel[:20])),'distinct_CDR3s':len({r['CDRs'][2] for r in panel}),'distinct_docking_pose_ids':len({r['pose_id'] for r in panel}),'maximum_entries_per_CDR3':max(collections.Counter(r['CDRs'][2] for r in panel).values()),'local_CDR3_maximum_edit_identity':max(r['local_novelty']['max_edit_identity'] for r in panel),'BLAST_databases_verified':list(dbs),'BLAST_final_query_coverage_complete':True,'ancestry_audit_all_pass':True,'format_roundtrip_pass':True,'selection_idempotent':True,'glycan_screen_explicit_for_all_final_and_potential_candidates':True,'all_three_biological_objectives':'UNVERIFIED','official_novelty_acceptance':'PENDING_ADAPTYV','batches':batches,'unbound_postprocessing_error_and_recovery':str(S/'reference/unbound_chain_audit_recovery.json'),'public_name_rng':'MT19937; seed 20061973; cosmetic naming only','files':frozen}
    (O/'release_validation.json').write_text(json.dumps(report,indent=2))
    (S/'reference/final_selection.json').write_text(json.dumps({'eligible_count':len(pool),'excluded_count':len(excluded),'selection_source':str(Path(__file__)),'selection_source_sha256':sha(__file__),'final_sequence_sha256s':[r['sequence_sha256'] for r in panel],'ranking':'Evidence stratum first, conservative pH proxy within stratum; 0.025 kcal per prior identical CDR3, plus 0.01 per prior pose for coarse-only reserves; maximum six per CDR3. Heuristic, not a probability or physical free-energy adjustment.'},indent=2))
    print(json.dumps({k:v for k,v in report.items() if k not in ['files','batches']},indent=2))
    for r in panel[:20]: print(r['rank'],r['public_name'],r['candidate_key'],r['tier'],round(r['ranking_contrast_kcal_proxy'],4))
if __name__=='__main__':main()
