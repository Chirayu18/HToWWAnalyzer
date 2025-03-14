import datetime
import vector
vector.register_awkward()
import os
import dask
import dask_awkward as dak
import awkward as ak
from coffea import processor
#from coffea.nanoevents.methods import candidate
from coffea.dataset_tools import (
    apply_to_fileset,
    preprocess,
)
from coffea.nanoevents import NanoAODSchema
from coffea.analysis_tools import PackedSelection
from argparse import ArgumentParser
from processing.LeptonSelections import *
parser = ArgumentParser()
parser.add_argument('--infile', default=None, help='/path/to/input/file')
parser.add_argument('--inputlist', default=None, help='/path/to/input/list')
parser.add_argument('--tag', default='', help='tag that is appended to output file name')
parser.add_argument('--outdir', default=None, help='/path/to/output/')
parser.add_argument('--islocal', default="True", help='boolean that decides whether to include xrootd redirector in dataset')
args = parser.parse_args()

cwd = os.getcwd()
if(args.inputlist != None):
    infile = open(args.inputlist, "r")
    paths = infile.read().splitlines()
    infile.close()
elif(args.infile != None):
     paths = [args.infile]

if not os.path.isdir(args.outdir): os.system("mkdir -p %s"%(args.outdir))

filesetname = args.tag

start_time = datetime.datetime.now()

print(start_time, ": Proceeding to run skim on", filesetname)

if(args.islocal == "True"):
    fileset = {
        filesetname: {
            "files": {"file://"+path: "Events" for path in paths},
        }
    }
else:
        fileset = {
        filesetname: {
            "files": {"root://cms-xrd-global.cern.ch/"+path: "Events" for path in paths},
        }
    }
dataset_runnable, dataset_updated = preprocess(
    fileset,
    align_clusters=False,
    step_size=100_000,
    skip_bad_files=True,
    save_form=False,
)

class GenInfo(processor.ProcessorABC):
    def __init__(self):
        pass

    def process(self, events):
        events, lead_trig_lep, sublead_trig_lep = applyTriggerPaths(events)
        print(lead_trig_lep.compute())
        print(sublead_trig_lep.compute())
        """
        genParticles = ak.zip(
        {
            "pt": events.GenPart.pt,
            "eta": events.GenPart.eta,
            "phi": events.GenPart.phi,
            "mass": events.GenPart.mass,
            "pdgId": events.GenPart.pdgId,
            "status": events.GenPart.status,
            "genPartIdxMother": events.GenPart.genPartIdxMother,
        },
        with_name="PtEtaPhiMCandidate",
        behavior=candidate.behavior,
        )
        """
        genParticles = ak.zip(
        {
            "pt": events.GenPart.pt,
            "eta": events.GenPart.eta,
            "phi": events.GenPart.phi,
            "mass": events.GenPart.mass,
            "pdgId": events.GenPart.pdgId,
            "status": events.GenPart.status,
            "genPartIdxMother": events.GenPart.genPartIdxMother,
        },
        with_name="LorentzVector"
        )

        muons = dak.zip(
        {
            "pt": events.Muon.pt,
            "eta": events.Muon.eta,
            "phi": events.Muon.phi,
            "mass": events.Muon.mass,
            "charge": events.Muon.charge,
            "isolation": events.Muon.pfRelIso04_all,
            "tightId": events.Muon.tightId,
            "mediumId": events.Muon.mediumId,
            "dxy": events.Muon.dxy,
            "dz": events.Muon.dz,
            "SIP3D": events.Muon.sip3d,
            "gen": events.Muon.genPartIdx
        },
        with_name="LorentzVector"
        )

        electrons = dak.zip(
        {
            "pt": events.Electron.pt,
            "eta": events.Electron.eta,
            "phi": events.Electron.phi,
            "mass": events.Electron.mass,
            "charge": events.Electron.charge,
            "cutBased": events.Electron.cutBased,
            "gen": events.Electron.genPartIdx
        },
        with_name="LorentzVector"
        )

        electrons = dak.zip(
        {
            "pt": events.Electron.pt,
            "eta": events.Electron.eta,
            "phi": events.Electron.phi,
            "mass": events.Electron.mass,
            "charge": events.Electron.charge,
            "cutBased": events.Electron.cutBased,
            "gen": events.Electron.genPartIdx
        },
        with_name="LorentzVector"
        )


        MET = dak.zip(
        {
            "pt": events.MET.pt,
            "phi": events.MET.phi,
        },
        with_name="LorentzVector"
        )

        electrons = electrons[dak.argsort(electrons.pt, ascending = False, axis=1)]
        #MET = MET[dak.argsort(MET.pt, ascending = False, axis=1)]
        muons = muons[dak.argsort(muons.pt, ascending = False, axis=1)]

        muons = muons[
            (muons.pt > 15)
            & (abs(muons.eta) < 2.4)
            & (muons.isolation < 0.15)
            #& (muons.SIP3D < 4)
            #& (muons.dz < 1)
            #& (muons.dxy < 0.5)
            & (muons.tightId == 1)
        ]

        electrons = electrons[
            (electrons.pt > 15)
            & (abs(electrons.eta) < 2.5)
            #& (electrons.isolation < 0.35)
            #& (electrons.SIP3D < 4)
            #& (electrons.dz < 1)
            #& (electrons.dxy < 0.5)
            & (electrons.cutBased >= 4)
        ]
    
        selection = PackedSelection()
        selection.add_multiple({
            "oneMuon": dak.num(muons) >= 1,
            "oneElectron": dak.num(electrons) >= 1,
            "minMET": MET.pt > 45,
        })    

        mask = selection.all("oneMuon","oneElectron","minMET")
        #mask = selection.all("minMET")


        genParticles["isHardProcess"] = events.GenPart.hasFlags('isHardProcess')
        genParticles["fromHardProcess"] = events.GenPart.hasFlags('fromHardProcess')
        genParticles["isLastCopy"] = events.GenPart.hasFlags('isLastCopy')
        genParticles["isLastCopyBeforeFSR"] = events.GenPart.hasFlags('isLastCopyBeforeFSR')
        #maskin = (genParticles["isHardProcess"]) & (genParticles["status"] == 21 )
        #maskfi = (genParticles["isHardProcess"]) & (genParticles["status"] > 21 )
        #initialState = genParticles[maskin]
        #finalState = genParticles[maskfi]
        #print(len(muons.compute()))
        muons = muons[mask]
        #print(len(muons.compute()))
        electrons = electrons[mask]
        genParticles = genParticles[mask]
        MET = MET[mask]
        out_dict = dak.zip({
                "genParticles": genParticles,
                "Muons": muons,
                "Electrons": electrons,
                "MET": MET,
        #        "initialState": initialState,
        #        "finalState": finalState,
        }, depth_limit=1)
        
        return dak.to_parquet(out_dict, destination=args.outdir, prefix = args.tag, compute=False) # Notice the compute=False flag

    def postprocess(self, accumulator):
        pass

to_compute = apply_to_fileset(
                GenInfo(),
                dataset_runnable,
                schemaclass=NanoAODSchema,
            )

(out, ) = dask.compute(to_compute)

end_time = datetime.datetime.now()

print(end_time, ": Sucessfully ran skim!")
