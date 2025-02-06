import os
import glob
import datetime
import uproot
import awkward as ak
import numpy as np
from argparse import ArgumentParser
import matplotlib.pyplot as plt

parser = ArgumentParser()
parser.add_argument('--inputlist', default=None, help='/path/to/input/dir')
parser.add_argument('--infile', default=None, help = 'path/to/input/file')
parser.add_argument('--tag', default=None, help='tag that determines output name. should correspond to process tag')
parser.add_argument('--outdir', default=None, help='/path/to/output/')
args = parser.parse_args()

if(args.inputlist != None):
    infile = open(args.inputlist, "r")
    files = infile.read().splitlines()
    infile.close()
if(args.infile != None):
    files = [args.infile]

def analyse(infile, outdir, tag, idx):

    events = ak.from_parquet(infile)

    genParticles = events.genParticles
    initialState = events.initialState
    finalState = events.finalState
    mask_Dc = ((abs(finalState.pdgId) > 10 ) & (abs(finalState.pdgId)<19)) | (abs(finalState.pdgId) == 24 )
    fs_noDc = finalState[~mask_Dc]
    inP = abs(initialState.pdgId)
    fP = abs(fs_noDc.pdgId)
    finalState_lep = finalState[ ( (abs(finalState.pdgId) == 11) | abs(finalState.pdgId == 13) ) ]
    finalState_nu = finalState[ ( (abs(finalState.pdgId) == 12) | abs(finalState.pdgId == 14) ) ]
    finalState_pairs = ak.cartesian([finalState_lep, finalState_nu])
    #finalState_pairs = finalState_pairs[ ak.num(finalState_pairs) > 0 ]
    finalState_mothIdx = ak.cartesian([finalState_lep.genPartIdxMother,finalState_nu.genPartIdxMother])
    w, i  = ak.unzip(finalState_mothIdx)
    condition = w==i
    finalState_pairs = finalState_pairs[condition]
    leps,nus = ak.unzip(finalState_pairs)
    leps = ak.flatten(leps)
    nus = ak.flatten(nus)
    delta_phi = leps.phi - nus.phi
    delta_phi = ak.where(delta_phi > np.pi, 2 * np.pi - delta_phi, delta_phi)
    mT = np.sqrt( 2 * leps.pt * nus.pt * ( 1 - np.cos(delta_phi)))
    plt.hist(mT, histtype='step',bins=70)
    plt.xlabel("Mt of W")
    plt.savefig("mT.png")
    plt.savefig("mT.pdf")
    print(mT)

    def checkProcess(initialState,finalState):
        valid_ini = ak.any(inP == initialState[0],axis=1) # Makes sure EVERY initial state particle is present at least once
        valid_ini2 = inP == initialState[0] # Makes sure every particle that exists in the initial state is present in our list of possible candidates
        if(len(initialState)>1):
            for i in initialState[1:]:
                valid_ini = (valid_ini) & (ak.any(inP == i,axis=1))
                valid_ini2 = (valid_ini2) | (inP ==i)
        valid_ini2 = ak.all(valid_ini2, axis=1)
        valid_ini = valid_ini & valid_ini2 & ( ak.num(inP) == len(initialState))
        valid_final = ak.any(fP == finalState[0],axis=1)
        valid_final2 = fP == finalState[0]
        if(len(finalState)>1):
            for i in finalState[1:]:
                valid_final = (valid_final) & (ak.any(fP == i,axis=1))
                valid_final2 = (valid_final2) | (fP ==i)
        valid_final2 = ak.all(valid_final2, axis=1)
        valid_final = valid_final & valid_final2 & ( ak.num(fP) == len(finalState))
        #valid_final = ak.all(valid_final, axis=1)
        return (valid_ini) & (valid_final)

    events_cch = events[checkProcess([ 4 , 4 ], [25])]
    print("Percentage of cc - > h: ", float(ak.num(events_cch,axis=0))/ak.num(events,axis=0))

    events_cchg = events[checkProcess([ 4 , 4 ], [25,21])]
    print("Percentage of cc - > hg: ", float(ak.num(events_cchg,axis=0))/ak.num(events,axis=0))

    events_cchgg = events[checkProcess([ 4 , 4 ], [25,21,21])]
    print("Percentage of cc - > hgg: ", float(ak.num(events_cchgg,axis=0))/ak.num(events,axis=0))

    events_gchc = events[checkProcess([ 21 , 4 ], [25,4])]
    print("Percentage of gc - > hc: ", float(ak.num(events_gchc,axis=0))/ak.num(events,axis=0))

    events_gghcc = events[checkProcess([ 21 , 21 ], [25,4,4])]
    print("Percentage of gg - > hcc: ", float(ak.num(events_gghcc,axis=0))/ak.num(events,axis=0))

    events_gcghc = events[checkProcess([21,4],[21,25,4]) | checkProcess([1,4],[1,25,4])| checkProcess([2,4],[2,25,4])| checkProcess([3,4],[3,25,4])]
    print("Percentage of l(g)c - > l(g)hc: ", float(ak.num(events_gcghc,axis=0))/ak.num(events,axis=0))

    outfile = outdir+"/"+tag+"-part"+str(idx)+".root"

    """
    with uproot.recreate(outfile) as fout:
            fout["Events"] = ({
                "muon_1": Higgs.muons.mu_1,
                "muon_2": Higgs.muons.mu_2,
                "muon_3": Higgs.muons.mu_3,
                "muon_4": Higgs.muons.mu_4,
                "Z1": Higgs.z1,
                "Z2": Higgs.z2,
                "Higgs": Higgs.Higgs,
                "jet": jet,
                "goodParton": goodParton,
                "goodJet": goodJet,
                "Weights": weights,
                "PSWeights": PSWeights,
                "LHEScaleWeights": LHEScaleWeights,
                "PHEPdfWeights": LHEPdfWeights
                })
    """
        
start_time = datetime.datetime.now()

print(start_time, ": Proceeding to run analysis on: \n", files)

for file, idx in zip(files, range(1, len(files)+1)):
    analyse(file, args.outdir, args.tag, idx)

end_time = datetime.datetime.now()

print(end_time, ": ran succesfully.")
