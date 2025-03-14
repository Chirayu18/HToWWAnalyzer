import getpass

filesPerJob = 2
processes = ["tqH_HToZZTo4L_M125_TuneCP5_13TeV-jhugenv7011-pythia8"]
tag = "skimtestd"
workdirBasePath = "/user/"+getpass.getuser()+"/NanoHJetAnalyser/HJetAnalyser/workdirs"
envPath = "/user/"+getpass.getuser()+"/NanoHJetAnalyser/HJetAnalyser/env"
outdir = "/pnfs/iihe/cms/store/user/"+getpass.getuser()+"/ntuples/NanoAOD/"
jobFlavour = "workday"
shift = "nominal"
era = "2018"


#The following list contains all backgrounds
#processes = ["bbH_HToZZTo4L_M125_TuneCP2_13TeV-jhugenv7011-pythia8", "DYJetsToLL_M-10to50_TuneCP5_13TeV-amcatnloFXFX-pythia8", "DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8", "GluGluHToZZTo4L_M125_TuneCP5_13TeV_powheg2_JHUGenV7011_pythia", "GluGluToContinToZZTo4e_TuneCP5_13TeV-mcfm701-pythia8", "GluGluToContinToZZTo4mu_TuneCP5_13TeV-mcfm701-pythia8", "GluGluToContinToZZTo4tau_TuneCP5_13TeV-mcfm701-pythia8", "GluGluToContinToZZTo2e2mu_TuneCP5_13TeV-mcfm701-pythia8", "GluGluToContinToZZTo2e2tau_TuneCP5_13TeV-mcfm701-pythia8", "GluGluToContinToZZTo2mu2tau_TuneCP5_13TeV-mcfm701-pythia8", "tqH_HToZZTo4L_M125_TuneCP5_13TeV-jhugenv7011-pythia8", "ttH_HToZZ_4LFilter_M125_TuneCP5_13TeV_powheg2_JHUGenV7011_pythia8", "VBF_HToZZTo4L_M125_TuneCP5_13TeV_powheg2_JHUGenV7011_pythia8", "WminusH_HToZZTo4L_M125_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV7011_pythia8", "WplusH_HToZZTo4L_M125_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV7011_pythia8","ZH_HToZZ_4LFilter_M125_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV7011_pythia8", "ZZTo4L_TuneCP5_13TeV_powheg_pythia8"]


"""
https://batchdocs.web.cern.ch/local/submit.html
espresso     = 20 minutes
microcentury = 1 hour
longlunch    = 2 hours
workday      = 8 hours
tomorrow     = 1 day
testmatch    = 3 days
nextweek     = 1 week
"""