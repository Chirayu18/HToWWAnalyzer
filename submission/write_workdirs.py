import submission_config as cfg
from argparse import ArgumentParser
import os
import json
import glob

parser = ArgumentParser()
parser.add_argument('--mode', default="", help='either "skim" or "analysis"')
parser.add_argument('--skimdir', default="", help='path/to/dir/of/skim. only relevant in analysis mode')
args = parser.parse_args()

if(args.mode == ""):
    raise Exception('Please specify a mode you would like to run!')

if os.path.isdir("%s/"%(cfg.workdirBasePath+"/"+cfg.tag+"/"+cfg.shift+"/"+args.mode)):
    print("these workdirs already existed. consider deleting them before running again for a clean setup...")
    print("would you like to delete them now? y/n")
    choice = input()
    if(choice == "y" or choice == "yes"):
        os.system("rm -rf "+cfg.workdirBasePath+"/"+cfg.tag+"/"+cfg.shift+"/"+args.mode)
        print("deleted existing workdirs...")
else:
    os.system("mkdir -p %s/"%(cfg.workdirBasePath+"/"+cfg.tag+"/"+cfg.shift+"/"+args.mode))

# Opening JSON file
jsonfile = open('HJetSamples/metadata/NanoAOD/2018_samples.json')
 
# returns JSON object as 
# a dictionary
sampleinfo = json.load(jsonfile)

if(args.mode == "skim"):
    for p in cfg.processes:
        filelist = "./HJetSamples/filelists/"+cfg.era+"/NanoAOD/"+p+".txt"
        print("Writing workdirs for files in: " + filelist)
        no_of_files = sum(1 for line in open(filelist))
        nlastjob = no_of_files % cfg.filesPerJob
        if(nlastjob == 0):
            njobs = int((no_of_files/cfg.filesPerJob))
        else: 
            njobs = int((no_of_files/cfg.filesPerJob) + 1)

        workdirs_flag = False
        #creates workdir for each job
        for i in range(1,njobs+1):
            if not os.path.isdir("%s/job_%s"%(cfg.workdirBasePath+"/"+cfg.tag+"/"+cfg.shift+"/"+args.mode+"/"+p, i)): os.system("mkdir -p %s/job_%s"%(cfg.workdirBasePath+"/"+cfg.tag+"/"+cfg.shift+"/"+args.mode+"/"+p, i))
            else: workdirs_flag = True

        #writes filelist from masterfile list for each job
        with open(filelist, 'r') as master:
            content = master.readlines()
            for i in range(1,njobs):
                f = open("%s/job_%s/file_list.txt"%(cfg.workdirBasePath+"/"+cfg.tag+"/"+cfg.shift+"/"+args.mode+"/"+p, i), 'w')
                for k in range(cfg.filesPerJob):
                    f.write(content[k + (i-1)*cfg.filesPerJob])
                f.close()

            f = open("%s/job_%s/file_list.txt"%(cfg.workdirBasePath+"/"+cfg.tag+"/"+cfg.shift+"/"+args.mode+"/"+p, njobs), 'w')

            if(nlastjob == 0):
                counter = cfg.filesPerJob
            else: 
                counter = no_of_files%cfg.filesPerJob
            for k in range(counter):
                f.write(content[k + (njobs-1)*cfg.filesPerJob])
            f.close()
            
        for i in range(1,njobs+1): 
            localdir_job = "%s/job_%s"%(cfg.workdirBasePath+"/"+cfg.tag+"/"+cfg.shift+"/"+args.mode+"/"+p, str(i))

            fskim_ = open(localdir_job+"/launch.sh", 'w')
            fskim_.write("#!/bin/bash \n")
            fskim_.write("pwd=$PWD \n")
            fskim_.write("cd /user/fheyen/NanoHJetAnalyser/HJetAnalyser \n")  
            fskim_.write("cd $pwd \n")
            fskim_.write("export X509_USER_PROXY=/user/$USER/x509up_u$(id -u $USER) \n")
            fskim_.write("export HOME=/user/$USER \n")
            fskim_.write("source %s/bin/activate \n"%(cfg.envPath))
            if(args.mode == "skim"):
                fskim_.write("python /user/fheyen/NanoHJetAnalyser/HJetAnalyser/skim/skim.py --inputlist %s --outdir %s --islocal %s --tag %s\n"%(localdir_job+"/file_list.txt", cfg.outdir+"/"+cfg.tag+"/"+cfg.shift+"/"+args.mode+"/"+p, sampleinfo[p]["local"], p+"-job_"+str(i)))
            if(args.mode == "analysis"):
                fskim_.write("python /user/fheyen/NanoHJetAnalyser/HJetAnalyser/analysis/analyse.py --inputlist %s --outdir %s --tag %s \n"%(localdir_job+"/file_list.txt", cfg.outdir+"/"+cfg.tag+"/"+cfg.shift+"/"+args.mode+"/"+p+"/job_"+str(i), p))
            fskim_.close()

            os.system("chmod u+x "+localdir_job+"/launch.sh")

        # write condor submission script
        localdir_process = (cfg.workdirBasePath+"/"+cfg.tag+"/"+cfg.shift+"/"+args.mode+"/"+p)
        fsubmit_ = open(localdir_process+"/job.submit", 'w')
        fsubmit_.write("Universe = vanilla \n")
        fsubmit_.write("Arguments = $(number) \n")
        fsubmit_.write("Executable = %s/job_$(number)/launch.sh \n"%(localdir_process))
        fsubmit_.write(" \n")
        fsubmit_.write("Error = %s/job_$(number)/job.err \n"%(localdir_process))
        fsubmit_.write("Output = %s/job_$(number)/job.out \n"%(localdir_process))
        fsubmit_.write("Log = %s/job_$(number)/job.log \n"%(localdir_process))
        fsubmit_.write(" \n")
        fsubmit_.write('+jobFlavour = "%s" \n' %cfg.jobFlavour)
        fsubmit_.write(" \n")
        fsubmit_.write('queue number from %s/fileNumbers.txt'%(localdir_process))
        fsubmit_.close()

        ffileNumbers = open(localdir_process+"/fileNumbers.txt", 'w')
        for n in range(1,njobs+1):
            ffileNumbers.write(str(n)+"\n")

        ffileNumbers.close()

        submissiondict = {
            "filesPerJob": cfg.filesPerJob,
            "processes": cfg.processes,
            "tag": cfg.tag,
            "workdirBasepath": cfg.workdirBasePath,
            "envPath": cfg.envPath,
            "outdir": cfg.outdir,
            "jobFlavour": cfg.jobFlavour,
            "shift": cfg.shift,
            "era": cfg.era,
        }

        fconfig_ = open(localdir_process+"/config.py", 'w')
        fconfig_.write("submissionDict = { \n")
        fconfig_.write("    'filesPerJob': %s,\n"%(cfg.filesPerJob))
        fconfig_.write("    'process': '%s',\n"%(p))
        fconfig_.write("    'tag': '%s',\n"%(cfg.tag))
        fconfig_.write("    'workdirPath': '%s',\n"%(cfg.workdirBasePath+"/"+cfg.tag+"/"+cfg.shift+"/"+p))
        fconfig_.write("    'envPath': '%s',\n"%(cfg.envPath))
        fconfig_.write("    'outdir': '%s',\n"%(cfg.outdir))
        fconfig_.write("    'jobFlavour': '%s',\n"%(cfg.jobFlavour))
        fconfig_.write("    'shift': '%s',\n"%(cfg.shift))
        fconfig_.write("    'era': '%s',\n"%(cfg.era))
        fconfig_.write("    'nJobs': %s,\n"%(njobs))
        fconfig_.write("    'mode': %s\n"%(args.mode))
        fconfig_.write("}")

        fconfig_.close()

    fprocesses_ = open(cfg.workdirBasePath+"/"+cfg.tag+"/"+cfg.shift+"/"+args.mode+"/"+"/processes.txt", 'w')
    for p in cfg.processes:
        fprocesses_.write(p+"\n")
    fprocesses_.close()
elif(args.mode == "analysis"):
    processListFile = open(args.skimdir+"/processes.txt", "r")
    processList = processListFile.read().splitlines()
    processListFile.close()

    for p in processList:
        processDir = args.skimdir+"/"+p
        processFileList = glob.glob(processDir+"/*.parquet")

        print("Writing workdirs for files in: ", processDir)
        no_of_files = len(processFileList)
        nlastjob = no_of_files % cfg.filesPerJob
        if(nlastjob == 0):
            njobs = int((no_of_files/cfg.filesPerJob))
        else: 
            njobs = int((no_of_files/cfg.filesPerJob) + 1)

        workdirs_flag = False
        #creates workdir for each job
        for i in range(1,njobs+1):
            if not os.path.isdir("%s/job_%s"%(cfg.workdirBasePath+"/"+cfg.tag+"/"+cfg.shift+"/"+args.mode+"/"+p, i)): os.system("mkdir -p %s/job_%s"%(cfg.workdirBasePath+"/"+cfg.tag+"/"+cfg.shift+"/"+args.mode+"/"+p, i))
            else: workdirs_flag = True

        #writes filelist from masterfile list for each job
        for i in range(1,njobs):
            f = open("%s/job_%s/file_list.txt"%(cfg.workdirBasePath+"/"+cfg.tag+"/"+cfg.shift+"/"+args.mode+"/"+p, i), 'w')
            for k in range(cfg.filesPerJob):
                f.write(processFileList[k + (i-1)*cfg.filesPerJob])
                f.write("\n")
            f.close()

        f = open("%s/job_%s/file_list.txt"%(cfg.workdirBasePath+"/"+cfg.tag+"/"+cfg.shift+"/"+args.mode+"/"+p, njobs), 'w')

        if(nlastjob == 0):
            counter = cfg.filesPerJob
        else: 
            counter = no_of_files%cfg.filesPerJob
        for k in range(counter):
            f.write(processFileList[k + (njobs-1)*cfg.filesPerJob])
            f.write("\n")
        f.close()
            
        for i in range(1,njobs+1): 
            localdir_job = "%s/job_%s"%(cfg.workdirBasePath+"/"+cfg.tag+"/"+cfg.shift+"/"+args.mode+"/"+p, str(i))

            fskim_ = open(localdir_job+"/launch.sh", 'w')
            fskim_.write("#!/bin/bash \n")
            fskim_.write("pwd=$PWD \n")
            fskim_.write("cd /user/fheyen/NanoHJetAnalyser/HJetAnalyser \n")  
            fskim_.write("cd $pwd \n")
            fskim_.write("export X509_USER_PROXY=/user/$USER/x509up_u$(id -u $USER) \n")
            fskim_.write("export HOME=/user/$USER \n")
            fskim_.write("source %s/bin/activate \n"%(cfg.envPath))
            if(args.mode == "skim"):
                fskim_.write("python /user/fheyen/NanoHJetAnalyser/HJetAnalyser/skim/skim.py --inputlist %s --outdir %s --islocal %s --tag %s\n"%(localdir_job+"/file_list.txt", cfg.outdir+"/"+cfg.tag+"/"+cfg.shift+"/"+args.mode+"/"+p, sampleinfo[p]["local"], p+"-job_"+str(i)))
            if(args.mode == "analysis"):
                fskim_.write("python /user/fheyen/NanoHJetAnalyser/HJetAnalyser/analysis/analyse.py --inputlist %s --outdir %s --tag %s \n"%(localdir_job+"/file_list.txt", cfg.outdir+"/"+cfg.tag+"/"+cfg.shift+"/"+args.mode+"/"+p+"/job_"+str(i), p))
            fskim_.close()

            os.system("chmod u+x "+localdir_job+"/launch.sh")

        # write condor submission script
        localdir_process = (cfg.workdirBasePath+"/"+cfg.tag+"/"+cfg.shift+"/"+args.mode+"/"+p)
        fsubmit_ = open(localdir_process+"/job.submit", 'w')
        fsubmit_.write("Universe = vanilla \n")
        fsubmit_.write("Arguments = $(number) \n")
        fsubmit_.write("Executable = %s/job_$(number)/launch.sh \n"%(localdir_process))
        fsubmit_.write(" \n")
        fsubmit_.write("Error = %s/job_$(number)/job.err \n"%(localdir_process))
        fsubmit_.write("Output = %s/job_$(number)/job.out \n"%(localdir_process))
        fsubmit_.write("Log = %s/job_$(number)/job.log \n"%(localdir_process))
        fsubmit_.write(" \n")
        fsubmit_.write('+jobFlavour = "%s" \n' %cfg.jobFlavour)
        fsubmit_.write(" \n")
        fsubmit_.write('queue number from %s/fileNumbers.txt'%(localdir_process))
        fsubmit_.close()

        ffileNumbers = open(localdir_process+"/fileNumbers.txt", 'w')
        for n in range(1,njobs+1):
            ffileNumbers.write(str(n)+"\n")

        ffileNumbers.close()

        submissiondict = {
            "filesPerJob": cfg.filesPerJob,
            "processes": cfg.processes,
            "tag": cfg.tag,
            "workdirBasepath": cfg.workdirBasePath,
            "envPath": cfg.envPath,
            "outdir": cfg.outdir,
            "jobFlavour": cfg.jobFlavour,
            "shift": cfg.shift,
            "era": cfg.era,
        }

        fconfig_ = open(localdir_process+"/config.py", 'w')
        fconfig_.write("submissionDict = { \n")
        fconfig_.write("    'filesPerJob': %s,\n"%(cfg.filesPerJob))
        fconfig_.write("    'process': '%s',\n"%(p))
        fconfig_.write("    'tag': '%s',\n"%(cfg.tag))
        fconfig_.write("    'workdirPath': '%s',\n"%(cfg.workdirBasePath+"/"+cfg.tag+"/"+cfg.shift+"/"+p))
        fconfig_.write("    'envPath': '%s',\n"%(cfg.envPath))
        fconfig_.write("    'outdir': '%s',\n"%(cfg.outdir))
        fconfig_.write("    'jobFlavour': '%s',\n"%(cfg.jobFlavour))
        fconfig_.write("    'shift': '%s',\n"%(cfg.shift))
        fconfig_.write("    'era': '%s',\n"%(cfg.era))
        fconfig_.write("    'nJobs': %s,\n"%(njobs))
        fconfig_.write("    'mode': %s\n"%(args.mode))
        fconfig_.write("}")

        fconfig_.close()

        fprocesses_ = open(cfg.workdirBasePath+"/"+cfg.tag+"/"+cfg.shift+"/"+args.mode+"/"+"/processes.txt", 'w')
        for p in cfg.processes:
            fprocesses_.write(p+"\n")
        fprocesses_.close()

if(args.mode == "skim"):
    if not(os.path.isdir("%s/"%(cfg.outdir+"/"+cfg.tag+"/"+cfg.shift+"/"+args.mode))):
        os.system("echo writing output dir at %s"%(cfg.outdir+"/"+cfg.tag+"/"+cfg.shift+"/"+args.mode))
        os.system("mkdir -p %s/"%(cfg.outdir+"/"+cfg.tag+"/"+cfg.shift+"/"+args.mode))
        os.system("cp %s %s"%(cfg.workdirBasePath+"/"+cfg.tag+"/"+cfg.shift+"/"+args.mode+"/"+"/processes.txt", cfg.outdir+"/"+cfg.tag+"/"+cfg.shift+"/"+args.mode+"/"+"/processes.txt"))