import os
import glob
import datetime
import uproot
import awkward as ak
import numpy as np
from argparse import ArgumentParser
import matplotlib.pyplot as plt
import vector
from scipy.stats import norm

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

def plotEff(data,xlabel,title,filename):
    
    # Plot histogram with normalization
    bins = 70
    counts, bin_edges, _ = plt.hist(data, bins=bins, density=True, alpha=0.6, color='blue', edgecolor="black", label="Data")

    # Fit a Gaussian (normal distribution) to the data
    mu, sigma = norm.fit(data)  # Get mean and std dev

    # Generate x values for the Gaussian curve
    x = np.linspace(bin_edges[0], bin_edges[-1], 1000)
    pdf = norm.pdf(x, mu, sigma)  # Compute Gaussian curve

    # Plot the fitted Gaussian curve
    plt.plot(x, pdf, 'r-', label=f'Gaussian Fit\n$\mu$ = {mu:.2f}, $\sigma$ = {sigma:.2f}')

    # Apply CMS-style formatting
    plt.xlabel(xlabel, fontsize=14)
    plt.ylabel("Normalized Entries", fontsize=14)
    plt.title(title, fontsize=16,)  # "CMS" Style title
    plt.grid(True, linestyle="--", alpha=0.5)
    plt.legend(fontsize=12, loc="upper left")
    plt.savefig(f"{filename}.png")
    plt.savefig(f"{filename}.pdf")

    #plt.show()

def analyse(infile, outdir, tag, idx):

    events = ak.from_parquet(infile)

    genParticles = events.genParticles
    electrons = events.Electrons
    muons = events.Muons
    met = events.MET

    gen_elec = genParticles[electrons.gen]
    gen_muon = genParticles[muons.gen]

    two_muons = ak.num(muons) > 1 
    two_elec = ak.num(electrons) > 1
    print("Parents of muons of events with multiple muons:")
    print(genParticles[two_muons][gen_muon[two_muons].genPartIdxMother].pdgId.show())
    genParticles = genParticles[~(two_muons | two_elec)]
    muons = muons[~(two_muons | two_elec)]
    electrons = electrons[~(two_muons | two_elec)]
    met = met[~(two_muons | two_elec)]
    gen_muon = gen_muon[~(two_muons | two_elec)]
    gen_elec = gen_elec[~(two_muons | two_elec)]

    getMother = lambda eve,x,pdg: genParticles[eve, x.genPartIdxMother ] if abs(genParticles[eve, x.genPartIdxMother].pdgId) != pdg else getMother(eve,genParticles[eve,x.genPartIdxMother],pdg)
    gen_elec_mothers = ak.Array([ getMother(eve,x,11) for eve,x in enumerate(gen_elec) ])#genParticles[ gen_elec.genPartIdxMother]
    gen_muon_mothers = ak.Array([ getMother(eve,x,13) for eve,x in enumerate(gen_muon) ])#genParticles[ gen_elec.genPartIdxMother]

    plt.hist(gen_elec_mothers.pdgId,bins=50,range=(-25,25))
    plt.title("PdgId of mother of matched gen electron(Before matching)")
    plt.savefig("plots/reco_elec_motherspdgId_nogenMatch.png")
    plt.savefig("plots/reco_elec_motherspdgId_nogenMatch.pdf")
    plt.clf()

    eventMask = ak.flatten((abs(gen_elec_mothers.pdgId) == 24) & (abs(gen_muon_mothers.pdgId) == 24))
    genParticles = genParticles[eventMask]
    electrons = electrons[eventMask]
    muons = muons[eventMask]
    met = met[eventMask]
    gen_elec_mothers = gen_elec_mothers[eventMask]
    gen_muon_mothers = gen_muon_mothers[eventMask]
    gen_elec = gen_elec[eventMask]
    gen_muon = gen_muon[eventMask]

    plt.hist(gen_elec_mothers.pdgId,bins=50,range=(-25,25))
    plt.title("PdgId of mother of matched gen electron(After matching)")
    plt.savefig("plots/reco_elec_motherspdgId.png")
    plt.savefig("plots/reco_elec_motherspdgId.pdf")
    plt.clf()

    print("Events: ",len(gen_elec))

    plt.hist(gen_elec.pt,histtype='step',range=(0,100),bins=80,label="Matched Gen Electron")
    plt.hist(electrons.pt,histtype='step',range=(0,100),bins=80,label="Reco Electron")
    plt.xlabel("pT")
    plt.title("Pt of Electrons")
    plt.legend()
    plt.savefig("plots/reco_elec_pt.png")
    plt.savefig("plots/reco_elec_pt.pdf")
    plt.clf()

    plotEff(np.clip(electrons.pt/gen_elec.pt - 1,-0.5,0.5),"(Reco_pt - Gen_pt)/Gen_pt","Reconstruction efficiency for electrons","plots/reco_eff_elec_pt")
    plt.clf()


    plt.hist(gen_elec.eta,histtype='step',bins=80,label="Matched Gen Electron")
    plt.hist(electrons.eta,histtype='step',bins=80,label="Reco Electron")
    plt.xlabel("eta")
    plt.title("Eta of Electrons")
    plt.legend()
    plt.savefig("plots/reco_elec_eta.png")
    plt.savefig("plots/reco_elec_eta.pdf")
    plt.clf()

    plotEff(muons.pt/gen_muon.pt - 1,"(Reco_pt - Gen_pt)/Gen_pt","Reconstruction efficiency for muons","plots/reco_eff_muon_pt")
    plt.clf()

    plt.hist(gen_muon.eta,histtype='step',bins=80,label="Matched Gen Muon")
    plt.hist(muons.eta,histtype='step',bins=80,label="Reco Muon")
    plt.title("Eta of Muons")
    plt.xlabel("eta")
    plt.legend()
    plt.savefig("plots/reco_muon_eta.png")
    plt.savefig("plots/reco_muon_eta.pdf")
    plt.clf()


    maskin = (genParticles.isHardProcess) & (genParticles.status == 21 )
    maskfi = (genParticles.isHardProcess) & (genParticles.status > 21 )
    initialState = genParticles[maskin]
    finalState = genParticles[maskfi]
    mask_Dc = ((abs(finalState.pdgId) > 10 ) & (abs(finalState.pdgId)<19)) | (abs(finalState.pdgId) == 24 )
    finalState = finalState[~mask_Dc]
    inP = abs(initialState.pdgId)
    fP = abs(finalState.pdgId)
    #leps = genParticles[  (abs(genParticles.pdgId) == 11) | (abs(genParticles.pdgId) == 13 ) ]
    ws = genParticles[ ( abs(genParticles.pdgId) == 24 ) & (genParticles.status==22) ]

    gen_muon = genParticles[muons.gen]
    gen_muon_mothers = genParticles[ gen_muon.genPartIdxMother]

    leps_rec = ak.concatenate([electrons,muons],axis=1)
    leps_rec = leps_rec[ak.argsort(leps_rec.pt,ascending=False)]
    leps_gen = genParticles[leps_rec.gen]

    lep1_rec_vec = vector.array({
        "pt": leps_rec[:,0].pt,
        "eta": leps_rec[:,0].eta,
        "phi": leps_rec[:,0].phi,
        "mass": leps_rec[:,0].mass
    })
    lep2_rec_vec = vector.array({
        "pt": leps_rec[:,1].pt,
        "eta": leps_rec[:,1].eta,
        "phi": leps_rec[:,1].phi,
        "mass": leps_rec[:,1].mass
    })
    lep1_gen_vec = vector.array({
        "pt": leps_gen[:,0].pt,
        "eta": leps_gen[:,0].eta,
        "phi": leps_gen[:,0].phi,
        "mass": leps_gen[:,0].mass
    })
    lep2_gen_vec = vector.array({
        "pt": leps_gen[:,1].pt,
        "eta": leps_gen[:,1].eta,
        "phi": leps_gen[:,1].phi,
        "mass": leps_gen[:,1].mass
    })

    m_ll_rec = (lep1_rec_vec + lep2_rec_vec).mass
    m_ll_gen = (lep1_gen_vec + lep2_gen_vec).mass
    plt.hist(m_ll_rec,histtype='step',bins=80,label="RECO")
    plt.axvline(x=12, color='r', linestyle='--', linewidth=2, label="Dilepton Selection(m_ll > 12 GeV)")
    plt.axvline(x=72, color='b', linestyle='--', linewidth=2, label="High mll Selection(m_ll < 72 GeV)")
    plt.hist(m_ll_gen,histtype='step',bins=80,label="Matched Gen")
    plt.title("m_ll")
    plt.xlabel("m_ll")
    plt.legend()
    plt.savefig("plots/reco_mll.png")
    plt.savefig("plots/reco_mll.pdf")
    plt.clf()

    pt_ll_rec = (lep1_rec_vec + lep2_rec_vec).pt
    pt_ll_gen = (lep1_gen_vec + lep2_gen_vec).pt
    plt.hist(pt_ll_rec,histtype='step',range=(0,150),bins=80,label="RECO")
    plt.axvline(x=30, color='r', linestyle='--', linewidth=2, label="Dilepton Selection(pt_ll > 30 GeV)")
    plt.hist(pt_ll_gen,histtype='step',range=(0,150),bins=80,label="Matched Gen")
    plt.title("pt_ll")
    plt.xlabel("pt_ll")
    plt.legend()
    plt.savefig("plots/reco_ptll.png")
    plt.savefig("plots/reco_ptll.pdf")
    plt.clf()

    delta_phi = lep1_rec_vec.phi - lep2_rec_vec.phi
    delta_phi = ak.where(delta_phi > np.pi, 2 * np.pi - delta_phi, delta_phi)
    dR_ll_rec = np.sqrt( (lep1_rec_vec.eta - lep2_rec_vec.eta)**2 + delta_phi**2 )
    delta_phi = lep1_gen_vec.phi - lep2_gen_vec.phi
    delta_phi = ak.where(delta_phi > np.pi, 2 * np.pi - delta_phi, delta_phi)
    dR_ll_gen = np.sqrt( (lep1_gen_vec.eta - lep2_gen_vec.eta)**2 + delta_phi**2 )
    plt.hist(dR_ll_rec,histtype='step',bins=80,label="RECO")
    plt.axvline(x=0.4, color='r', linestyle='--', linewidth=2, label="Dilepton Selection(dR > 0.4)")
    plt.hist(dR_ll_gen,histtype='step',bins=80,label="Matched Gen")
    plt.title("dR(l1,l2)")
    plt.xlabel("dR")
    plt.legend()
    plt.savefig("plots/reco_dRll.png")
    plt.savefig("plots/reco_dRll.pdf")
    plt.clf()

    delta_phi = (lep1_rec_vec + lep2_rec_vec).phi - met.phi
    delta_phi = ak.where(delta_phi > np.pi, 2 * np.pi - delta_phi, delta_phi)
    mT_ll_rec = np.sqrt( 2 * pt_ll_rec * met.pt * ( 1 - np.cos(delta_phi)))
    delta_phi = (lep1_gen_vec + lep2_gen_vec).phi - met.phi
    delta_phi = ak.where(delta_phi > np.pi, 2 * np.pi - delta_phi, delta_phi)
    mT_ll_gen = np.sqrt( 2 * pt_ll_gen * met.pt * ( 1 - np.cos(delta_phi)))
    plt.hist(mT_ll_rec,histtype='step',bins=80,label="RECO")
    plt.axvline(x=60, color='r', linestyle='--', linewidth=2, label="mT Selection(mT > 60 GeV)")
    plt.hist(mT_ll_gen,histtype='step',bins=80,label="Matched Gen")
    plt.title("mT of dileptons")
    plt.xlabel("mT(ll)")
    plt.legend()
    plt.savefig("plots/reco_mTll.png")
    plt.savefig("plots/reco_mTll.pdf")
    plt.clf()
    #print(mT_ll_rec)

    delta_phi = (lep1_rec_vec).phi - met.phi
    delta_phi = ak.where(delta_phi > np.pi, 2 * np.pi - delta_phi, delta_phi)
    mT_l1_rec = np.sqrt( 2 * lep1_rec_vec.pt * met.pt * ( 1 - np.cos(delta_phi)))

    delta_phi = (lep2_rec_vec).phi - met.phi
    delta_phi = ak.where(delta_phi > np.pi, 2 * np.pi - delta_phi, delta_phi)
    mT_l2_rec = np.sqrt( 2 * lep2_rec_vec.pt * met.pt * ( 1 - np.cos(delta_phi)))
    delta_phi = (lep2_gen_vec).phi - met.phi
    delta_phi = ak.where(delta_phi > np.pi, 2 * np.pi - delta_phi, delta_phi)
    mT_l2_gen = np.sqrt( 2 * lep2_gen_vec.pt * met.pt * ( 1 - np.cos(delta_phi)))
    plt.hist(mT_l2_rec,histtype='step',bins=80,label="RECO")
    plt.axvline(x=30, color='r', linestyle='--', linewidth=2, label="mT Selection(mT > 30 GeV)")
    plt.hist(mT_l2_gen,histtype='step',bins=80,label="Matched Gen")
    plt.title("mT of sublead lepton")
    plt.xlabel("mT(l2)")
    plt.legend()
    plt.savefig("plots/reco_mTl2.png")
    plt.savefig("plots/reco_mTl2.pdf")
    plt.clf()

    plt.hist(leps_rec[:,0].charge*leps_rec[:,1].charge,histtype='step',range=(-1.5,-0.5))
    plt.title("Charge of lep1*lep2")
    plt.xlabel("charge*charge")
    plt.savefig("plots/reco_charge.png")
    plt.savefig("plots/reco_charge.pdf")
    plt.clf()
    
    """
    print(ws.mass)
    w1,w2 = ws[:,0], ws[:,1]
    w1_vec = vector.array({
    "pt": w1.pt,
    "eta": w1.eta,
    "phi": w1.phi,
    "mass": w1.mass
    })
    w2_vec = vector.array({
    "pt": w2.pt,
    "eta": w2.eta,
    "phi": w2.phi,
    "mass": w2.mass
    })
    higgs = w1_vec+w2_vec
    #print(w1_vec.mass)
    #print(higgs.mass)
    """

    #nus = genParticles[ ( (abs(finalState.pdgId) == 12) | abs(finalState.pdgId == 14) ) ]
    #lep_mot = [ [genParticles[event,j].mass for j in i] for event,i in enumerate(leps.genPartIdxMother) ]
    #print(lep_mot)
    #print(ws.mass)
    finalState_lep = genParticles[(abs(genParticles.pdgId) == 11) | (abs(genParticles.pdgId) == 13 )| (abs(genParticles.pdgId) == 15 ) ]
    finalState_nu = genParticles[(abs(genParticles.pdgId) == 12) | (abs(genParticles.pdgId) == 14  ) | (abs(genParticles.pdgId) == 16 )  ]
    leps_mothers = genParticles[ finalState_lep.genPartIdxMother ]
    nus_mothers = genParticles[ finalState_nu.genPartIdxMother ]

    leps = finalState_lep[ abs(leps_mothers.pdgId) == 24 ]
    leps_mothers = leps_mothers[ abs(leps_mothers.pdgId) == 24 ]
    sort = ak.argsort(leps_mothers.mass)
    leps_mothers = leps_mothers[sort]
    leps = leps[sort]

    nus = finalState_nu[ abs(nus_mothers.pdgId) == 24 ]
    nus_mothers = nus_mothers[ abs(nus_mothers.pdgId) == 24 ]
    sort = ak.argsort(nus_mothers.mass)
    nus_mothers = nus_mothers[sort]
    nus = nus[sort]
    print(leps.pdgId.show())

    lep1 = leps[:,0]
    lep1_vec = vector.array({
    "pt": lep1.pt,
    "eta": lep1.eta,
    "phi": lep1.phi,
    "mass": lep1.mass
    })
    nu1 = nus[:,0]
    nu1_vec = vector.array({
    "pt": nu1.pt,
    "eta": nu1.eta,
    "phi": nu1.phi,
    "mass": nu1.mass
    })
    w1 = lep1_vec + nu1_vec

    lep2 = leps[:,1]
    lep2_vec = vector.array({
    "pt": lep2.pt,
    "eta": lep2.eta,
    "phi": lep2.phi,
    "mass": lep2.mass
    })
    nu2 = nus[:,1]
    nu2_vec = vector.array({
    "pt": nu2.pt,
    "eta": nu2.eta,
    "phi": nu2.phi,
    "mass": nu2.mass
    })
    w2 = lep2_vec + nu2_vec
    delta_phi = lep1_vec.phi - nu1_vec.phi
    delta_phi = ak.where(delta_phi > np.pi, 2 * np.pi - delta_phi, delta_phi)
    mT1 = np.sqrt( 2 * lep1_vec.pt * nu1_vec.pt * ( 1 - np.cos(delta_phi)))
    delta_phi = lep2_vec.phi - nu2_vec.phi
    delta_phi = ak.where(delta_phi > np.pi, 2 * np.pi - delta_phi, delta_phi)
    mT2 = np.sqrt( 2 * lep2_vec.pt * nu2_vec.pt * ( 1 - np.cos(delta_phi)))

    plt.hist2d(lep1_vec.pt,lep1_vec.eta,range=[[0,150],[-6,6]],bins=50,cmin=1)
    plt.xlabel("Pt")
    plt.ylabel("Eta")
    plt.title("Lepton from Off Shell W")
    plt.savefig("plots/pt_eta_lep_offShell.png")
    plt.savefig("plots/pt_eta_lep_offShell.pdf")
    plt.clf()

    plt.hist2d(lep2_vec.pt,lep2_vec.eta,range=[[0,150],[-6,6]],bins=50,cmin=1)
    plt.xlabel("Pt")
    plt.ylabel("Eta")
    plt.title("Lepton from On Shell W")
    plt.savefig("plots/pt_eta_lep_onShell.png")
    plt.savefig("plots/pt_eta_lep_onShell.pdf")
    plt.clf()

    met_gen = nu1_vec + nu2_vec
    plt.hist2d(met_gen.pt,met_gen.eta,range=[[0,150],[-6,6]],bins=50,cmin=1)
    plt.xlabel("Pt")
    plt.ylabel("Eta")
    plt.title("Neutrino 1 + Neutrino 2")
    plt.savefig("plots/pt_eta_nu.png")
    plt.savefig("plots/pt_eta_nu.pdf")
    plt.clf()

    plt.hist(leps[:,0].pdgId,bins=60,histtype='step',label="Lepton from on Shell W")
    plt.hist(leps[:,1].pdgId,bins=60,histtype='step',label="Lepton from off Shell W")
    plt.xlabel("PdgID")
    plt.legend()
    plt.savefig("plots/pdgID_lep.png")
    plt.savefig("plots/pdgID_lep.pdf")
    plt.clf()

    plt.hist(nus[:,0].pdgId,bins=60,histtype='step',label="Neutrino from on Shell W")
    plt.hist(nus[:,1].pdgId,bins=60,histtype='step',label="Neutrino from off Shell W")
    plt.xlabel("PdgID")
    plt.legend()
    plt.savefig("plots/pdgID_nu.png")
    plt.savefig("plots/pdgID_nu.pdf")
    plt.clf()

    plt.hist(w1.mass, histtype='step',bins=70,label="Off Shell W")
    plt.hist(w2.mass, histtype='step',bins=70,label="On Shell W")
    plt.xlabel("Invariant mass of W")
    plt.legend()
    plt.savefig("plots/mass_W.png")
    plt.savefig("plots/mass_W.pdf")
    plt.clf()

    plt.hist(mT1, histtype='step',bins=70,label="Off Shell W")
    plt.hist(mT2, histtype='step',bins=70,label="On Shell W")
    plt.xlabel("Transverse mass (mT) of W")
    plt.legend()
    plt.savefig("plots/mT_W.png")
    plt.savefig("plots/mT_W.pdf")
    plt.clf()

    plt.hist((w1+w2).mass, histtype='step',bins=70,range=[120,130])
    plt.xlabel("Invariant mass of Higgs")
    plt.savefig("plots/mass_H.png")
    plt.savefig("plots/mass_H.pdf")
    plt.clf()

    """
    print(mT)
    #print(leps_mothers_onShell.pt*np.cosh(leps_mothers_onShell.eta).show())
    exit
    finalState_pairs = ak.cartesian([finalState_lep, finalState_nu])
    #finalState_pairs = finalState_pairs[ ak.num(finalState_pairs) > 0 ]
    finalState_mothIdx = ak.cartesian([finalState_lep.genPartIdxMother,finalState_nu.genPartIdxMother])
    w, i  = ak.unzip(finalState_mothIdx)
    condition = w==i
    w = w[condition]
    i = i[condition]
    finalState_pairs = finalState_pairs[condition]
    #print(leps_mothers.pdgId)
    mask_wmother = abs(leps_mothers.pdgId ) == 24 
    mask_onShell = leps_mothers.mass > 80
    finalState_pairs = finalState_pairs[mask_onShell]
    leps,nus = ak.unzip(finalState_pairs)
    leps = ak.flatten(leps)
    nus = ak.flatten(nus)
    delta_phi = leps.phi - nus.phi
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
    mT = np.sqrt( 2 * leps.pt * nus.pt * ( 1 - np.cos(delta_phi)))
    #print(mT)
    plt.hist(mT, histtype='step',bins=70)
    plt.xlabel("Mt of W")
    plt.savefig("mT.png")
    plt.savefig("mT.pdf")
    """

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
