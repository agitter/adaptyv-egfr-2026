"""Cutoff-sensitivity audit with the same receptor His set in related designs."""
from proton_ensemble import *

def main():
    out=S3/'intermediate/common_sites';out.mkdir(exist_ok=True)
    for cid in ['A00002','A00008']:
        probe(S3/'intermediate/refined'/(cid+'_human6ARU_relaxed.pdb'),out/(cid+'_human6ARU_common'),('A',),{i+1:i+334 for i in range(172)},[('A','37'),('A','50'),('A','100')])
    # Earlier parent had two adverse binder histidines in addition to the
    # three receptor sites. This preserves all sites, rather than truncating.
    probe(R/'stage2/intermediate/mm_binding/P00024_07_human6ARU_complex_relaxed.pdb',out/'P00024_07_human6ARU_common',('A',),{i+1:i+334 for i in range(172)},[('A','37'),('A','50'),('A','100'),('B','108'),('B','109')])
if __name__=='__main__':main()
