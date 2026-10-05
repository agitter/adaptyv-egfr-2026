"""Include the four tested domain-III receptor histidines for both lead designs.
This is a site-cutoff sensitivity test, not a pKa prediction.
"""
from proton_ensemble import *

def main():
    out=S3/'intermediate/common_sites';out.mkdir(exist_ok=True)
    for cid in ['C00003','C00009']:
        probe(S3/'intermediate/refined'/(cid+'_human6ARU_relaxed.pdb'),out/(cid+'_human6ARU_four_target_sites'),('A',),{i+1:i+334 for i in range(172)},[('A','25'),('A','37'),('A','50'),('A','100')])
if __name__=='__main__':main()
