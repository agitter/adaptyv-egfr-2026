"""Directional two-nitrogen histidine-clamp audit (classical geometry only)."""
from pathlib import Path
import sys,json,math
import numpy as np
R=Path('/mnt/data/egfr_campaign');S3=R/'stage3'
sys.path[:0]=[str(R/'stage2/code'),str(R/'code')]
from classical_design import unit,CHARGE_NAMES
TARGET=json.loads((R/'intermediate/human6ARU_aligned.json').read_text())
HIST=[r for r in TARGET if r['aa']=='H' and 334<=r['human_pos']<=505]

def option_features(aa,atoms):
    features=np.zeros((len(HIST),2))
    if aa not in 'DE':return features
    c=np.asarray(atoms['CG' if aa=='D' else 'CD'])
    for hi,h in enumerate(HIST):
        ha=h['atoms']
        for ni,(n,neighbors) in enumerate([('ND1',['CG','CE1']),('NE2',['CD2','CE1'])]):
            p=np.asarray(ha[n]);direction=unit(p-np.mean([ha[a] for a in neighbors],axis=0))
            for oxygen in CHARGE_NAMES[aa]:
                q=np.asarray(atoms[oxygen]);delta=q-p;dist=np.linalg.norm(delta)
                if not 2.3<dist<4.:continue
                donor=np.dot(direction,delta)/dist;acceptor=np.dot(unit(q-c),-delta)/dist
                value=math.exp(-((dist-2.85)/.45)**2)*np.clip((donor-.2)/.65,0,1)*np.clip((acceptor-.1)/.65,0,1)
                features[hi,ni]=max(features[hi,ni],value)
    return features

def audit(rows):
    feats=np.array([option_features(r['aa'],r['atoms']) for r in rows]);tot=feats.max(0);out=[]
    for h,f in zip(HIST,tot):
        hi=HIST.index(h)
        positions=[[rows[i].get('position',i+1) for i in np.where(feats[:,hi,n]>.2)[0]] for n in range(2)]
        out.append({'target_histidine':h['human_pos'],'ND1_best':float(f[0]),'NE2_best':float(f[1]),'two_sided_minimum':float(min(f)),'binder_positions_by_N':positions})
    return {'histidines':out,'strong_clamps':sum(h['two_sided_minimum']>.3 for h in out),'clamp_score':sum(h['two_sided_minimum'] for h in out)}

if __name__=='__main__':
    result=[]
    for p in (S3/'intermediate/designs').glob('*.json'):
        if p.stem.endswith('_summary'):continue
        d=json.loads(p.read_text());a=audit(d['structure']);result.append({'candidate_id':d['candidate_id'],'branch':d['branch'],**a})
    result.sort(key=lambda a:a['clamp_score'],reverse=True);(S3/'reference/clamp_audit.json').write_text(json.dumps(result,indent=2))
    for a in result[:12]:print(a['candidate_id'],'clamps',a['strong_clamps'],'score',round(a['clamp_score'],3),[h for h in a['histidines'] if h['ND1_best']+h['NE2_best']>.1])
