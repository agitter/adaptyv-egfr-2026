"""Bounded concurrent refinement; explicit per-candidate status and logs."""
from pathlib import Path
import sys,json,subprocess,time,os
R=Path('/mnt/data/egfr_campaign');S=R/'stage3'
selection=json.loads(Path(sys.argv[1]).read_text());workers=int(sys.argv[2]) if len(sys.argv)>2 else 2
active=[];done=[];pending=list(selection);start=time.time()
while active or pending:
    while pending and len(active)<workers:
        item=pending.pop(0)
        if isinstance(item,str):candidate=item;species='human6ARU'
        else:candidate=item['candidate_id'];species=item.get('species','human6ARU')
        log=S/'logs'/('mm_'+candidate+'_'+species+'.log');f=log.open('w')
        p=subprocess.Popen([sys.executable,str(S/'code/refine_new.py'),candidate,species],stdout=f,stderr=subprocess.STDOUT,env=os.environ.copy())
        active.append((p,f,candidate,species,str(log),time.time()));print('started',candidate,species,'pid',p.pid,flush=True)
    for item in active[:]:
        p,f,candidate,species,log,when=item;rc=p.poll()
        if rc is None:continue
        f.close();active.remove(item);record={'candidate_id':candidate,'species':species,'returncode':rc,'log':log,'seconds':time.time()-when};done.append(record);print('finished',candidate,species,rc,round(record['seconds'],1),flush=True)
        (S/'reference'/(Path(sys.argv[1]).stem.replace('_selection','')+'_status.json')).write_text(json.dumps({'completed':done,'pending':pending,'active':[{'candidate_id':a[2],'species':a[3],'pid':a[0].pid} for a in active],'seconds':time.time()-start},indent=2))
    time.sleep(1)
print('ALL FINISHED',len(done),'seconds',round(time.time()-start,1),flush=True)
