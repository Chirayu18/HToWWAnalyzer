import os
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('--indir', default=None, help='/path/to/dir')
args = parser.parse_args()

with open(args.indir+"/processes.txt") as file: 
    processes = file.read().splitlines()

for p in processes:
    print("-----------------------")
    print("Submitting jobs for "+p)
    print("-----------------------")

    os.system("condor_submit %s/job.submit"%(args.indir+"/"+p))