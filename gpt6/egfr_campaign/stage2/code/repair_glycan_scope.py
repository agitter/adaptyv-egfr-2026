"""Assign sugars by covalent connectivity, not by an assumed chain label."""
from pathlib import Path
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
import json,collections
R=Path('/mnt/data/egfr_campaign');S=R/'stage2'
c=MMCIF2Dict(str(R/'inputs/provided/6ARU.cif'));g=json.loads((R/'intermediate/human6ARU_glycans.json').read_text());sugars={(x['chain'],str(x['resid'])) for x in g};edges=collections.defaultdict(set);roots=set()
for i,kind in enumerate(c['_struct_conn.conn_type_id']):
 if kind!='covale':continue
 nodes=[(c[f'_struct_conn.ptnr{k}_auth_asym_id'][i],c[f'_struct_conn.ptnr{k}_auth_seq_id'][i]) for k in [1,2]]
 a,b=nodes
 if a in sugars and b in sugars:edges[a].add(b);edges[b].add(a)
 elif a[0]=='A' and a not in sugars and b in sugars:roots.add(b)
 elif b[0]=='A' and b not in sugars and a in sugars:roots.add(a)
seen=set(roots);front=list(roots)
while front:
 a=front.pop()
 for b in edges[a]-seen:seen.add(b);front.append(b)
kept=[a for a in g if (a['chain'],str(a['resid'])) in seen];excluded=[a for a in g if (a['chain'],str(a['resid'])) not in seen]
(S/'intermediate/target_glycans.json').write_text(json.dumps(kept));(S/'reference/glycan_scope_audit.json').write_text(json.dumps({'roots':sorted(roots),'target_sugar_residues':sorted(seen),'retained_atoms':len(kept),'excluded_atoms':len(excluded),'excluded_residues':sorted(set((a['chain'],a['resid']) for a in excluded)),'reason':'Only sugars covalently linked to EGFR chain A are target obstacles; antibody-attached glycan excluded.'},indent=2))
print('Target glycan atoms',len(kept),'other glycan atoms',len(excluded))
