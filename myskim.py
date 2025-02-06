import datetime
import os
import dask
import dask_awkward as dak
import awkward as ak
from coffea import processor
from coffea.nanoevents.methods import candidate
from coffea.dataset_tools import (
    apply_to_fileset,
    preprocess,
)
from coffea.nanoevents import NanoAODSchema
from coffea.analysis_tools import PackedSelection
from argparse import ArgumentParser
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
    
        genParticles["isHardProcess"] = events.GenPart.hasFlags('isHardProcess')
        genParticles["fromHardProcess"] = events.GenPart.hasFlags('fromHardProcess')
        genParticles["isLastCopy"] = events.GenPart.hasFlags('isLastCopy')
        genParticles["isLastCopyBeforeFSR"] = events.GenPart.hasFlags('isLastCopyBeforeFSR')
        maskin = (genParticles["isHardProcess"]) & (genParticles["status"] == 21 )
        maskfi = (genParticles["isHardProcess"]) & (genParticles["status"] > 21 )
        initialState = genParticles[maskin]
        finalState = genParticles[maskfi]
        out_dict = dak.zip({
                "genParticles": genParticles,
                "initialState": initialState,
                "finalState": finalState,
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
