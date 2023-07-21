import argparse
import time
import numpy as np
import ROOT as r
import tqdm


SIGNAL_PDG_ID = 1000006
MAX_ETA = 2.3
N_JET_MAX = 12
N_FEAT = 14
N_PART_PER_JET = 10
DELTA_R_MATCH = 0.4

r.gROOT.SetBatch(1)

## Unicode Shortcuts
# phi = \u03A6
# eta = \u03B7
# delta = \u0394

## Effic Functions
# Jet Triggers
def singleJetEffic(inputlist, ptVal, etaVal):
     num = 0
     ver = inputlist
     for i in range(len(ver)):
          if len(ver[i]) > 0:
               for j in range(len(ver[i])):
                    if (ver[i][j].Pt() > ptVal) and (abs(ver[i][j].Eta()) < etaVal):
                         num += 1
                         break
     print(f'The effic of single jets with pT > {ptVal} is {num / len(ver)} over {len(ver)} events')

def doubleJetEffic(inputlist, ptVal, etaVal):
     num = 0
     ver = inputlist
     for i in range(len(ver)):
          if len(ver[i]) >= 2:
               for j in range(len(ver[i])):
                    if (ver[i][j].Pt() > ptVal) and (abs(ver[i][j].Eta()) < etaVal):
                         for k in range(len(ver[i])):
                              if (k!=j) and (ver[i][k].Pt() > ptVal):
                                   num += 1
                                   break
                         break     
     print(f'The effic of double jets with pT > {ptVal} and |\u03B7| < {etaVal} is {num / len(inputlist)} over {len(inputlist)} events')

def doubleJetdeltaEtaEffic(inputlist, ptVal, etaVal, deltaEta):
     num = 0
     ver = inputlist
     for i in range(len(ver)):
          if len(ver[i]) >= 2:
               for j in range(len(ver[i])):
                    if (ver[i][j].Pt() > ptVal) and (abs(ver[i][j].Eta()) < etaVal):
                         for k in range(len(ver[i])):
                              if (k!=j) and (ver[i][k].Pt() > ptVal) and (abs(ver[i][k].Eta()) < etaVal) and (abs(abs(ver[i][j].Eta()) - abs(ver[i][k].Eta())) < deltaEta):
                                   num += 1
                                   break
                         break
     print(f'The effic of double jets + \u0394\u03B7 with pT > {ptVal} and |\u03B7| < {etaVal} and \u0394\u03B7 < {deltaEta} is {num / len(inputlist)} over {len(inputlist)} events')

def doubleJetMassEffic(inputlist, ptVal1, ptVal2, massVal, etaVal):
     num = 0
     ver = inputlist
     for i in range(len(ver)):
          for j in range(len(ver[i])):
               if (ver[i][j].Pt() >= ptVal1) and (abs(ver[i][j].Eta()) < etaVal):
                    for k in range(len(ver[i])):
                         if (ver[i][k].Pt() >= ptVal2) and ((ver[i][j].M() + ver[i][k].M()) > massVal) and (abs(ver[i][j].Eta()) < etaVal) and (k!=j):
                              num += 1
                              break
                    break
     print(f'The effic of double jets + mass with pT > {ptVal1}, {ptVal2}; two jets pT > {ptVal2} and M_jj > {massVal} is {num / len(inputlist)} over {len(inputlist)} events')

def doubleJetMass2Effic(inputlist, ptVal, etaVal, deltaEta, massVal):
     num = 0
     ver = inputlist
     for i in range(len(ver)):
          for j in range(len(ver[i])):
               if (ver[i][j].Pt() >= ptVal) and (abs(ver[i][j].Eta()) < etaVal):
                    for k in range(len(ver[i])):
                         if (ver[i][k].Pt() >= ptVal) and (abs(ver[i][k].Eta()) < etaVal) and (abs(abs(ver[i][j].Eta()) - abs(ver[i][k].Eta())) < deltaEta) and ((ver[i][j].M() + ver[i][k].M()) > massVal) and (k!=j):
                              num += 1
                              break
                    break
     print(f'The effic of double jets + mass with pT > {ptVal} and |\u03B7| < {etaVal} and \u0394\u03B7 < {deltaEta} and M_jj > {massVal} is {num / len(inputlist)} over {len(inputlist)} events')

def tripleJetEffic(inputlist, ptVal1, ptVal2, ptVal3, etaVal):
     num = 0
     ver = inputlist
     for i in range(len(ver)):
          for j in range(len(ver[i])):
               if ver[i][j].Pt() >= ptVal1:
                    for k in range(len(ver[i])):
                         if (ver[i][k].Pt() >= ptVal2) and (abs(ver[i][k].Eta()) < etaVal) and (k!=j):
                              for m in range(len(ver[i])):
                                   if (ver[i][m].Pt() >= ptVal3) and (abs(ver[i][k].Eta()) < etaVal) and (m!=k) and (m!=j):
                                        num += 1
                                        break
                              break
                    break
     print(f'The effic of triple jets with pT > {ptVal1}, {ptVal2}, {ptVal3}; two jets pT > {ptVal2}, {ptVal3} and |\u03B7| < {etaVal} is {num / len(inputlist)} over {len(inputlist)} events')         

# Energy Sum Triggers
def EtmissEffic(inputlist, EVal, etaVal):
     num = 0
     ver = inputlist
     for i in range(len(ver)):
          vectorsumX = np.zeros(1, dtype=object)
          vectorsumY = np.zeros(1, dtype=object)
          for j in range(len(ver[i])):
               if (abs(ver[i][j].Eta()) < etaVal):
                    vectorsumX[0] = vectorsumX[0] + ver[i][j].Px()
                    vectorsumY[0] = vectorsumX[0] + ver[i][j].Px()
          vectorsum = np.sqrt(vectorsumX[0]**2 + vectorsumY[0]**2)
         # print(f'Vector sum = {vectorsum} for event {i}')
          if (vectorsum > EVal):
               num += 1
     print(f'The effic of E_miss_t energy sum > {EVal} of jets with |\u03B7| < {etaVal} is {num / len(inputlist)} over {len(inputlist)} events')

def HtEffic(inputlist, EVal, ptVal, etaVal):
     num = 0
     ver = inputlist
     for i in range(len(ver)):
          scalarsum = np.zeros(1, dtype=object)
          for j in range(len(ver[i])):
               if (ver[i][j].Pt() > ptVal) and (abs(ver[i][j].Eta()) < etaVal):
                    scalarsum[0] = scalarsum[0] + ver[i][j].Pt()
         # print(f'Scalar sum = {scalarsum} for event {i}')
          if (scalarsum[0] > EVal):
               num += 1
     print(f'The effic of H_t energy sum > {EVal} of jets with pT > {ptVal} and |\u03B7| < {etaVal} is {num / len(inputlist)} over {len(inputlist)} events')

def EtEffic(inputlist, EVal, etaVal):
     num = 0
     ver = inputlist
     for i in range(len(ver)):
          scalarsum = np.zeros(1, dtype=object)
          for j in range(len(ver[i])):
               if (abs(ver[i][j].Eta()) < etaVal):
                    scalarsum[0] = scalarsum[0] + ver[i][j].Pt()
          if (scalarsum[0] > EVal):
               num += 1
     print(f'The effic of E_t energy sum > {EVal} of jets with |\u03B7| < {etaVal} is {num / len(inputlist)} over {len(inputlist)} events')

## Effic Plotters
def singleJetEfficVal(inputlist, ptVal):
     num = 0
     ver = inputlist
     for i in range(len(ver)):
          if len(ver[i]) > 0:
               for j in range(len(ver[i])):
                    if (ver[i][j].Pt() > ptVal) and (abs(ver[i][j].Eta()) < 2.4):
                         num += 1
                         break
     return (num / len(ver))
def doubleJetEfficVal(inputlist, ptVal):
     num = 0
     ver = inputlist
     for i in range(len(ver)):
          if len(ver[i]) >= 2:
               for j in range(len(ver[i])):
                    if (ver[i][j].Pt() > ptVal) and (abs(ver[i][j].Eta()) < 2.4):
                         for k in range(len(ver[i])):
                              if (k!=j) and (ver[i][k].Pt() > ptVal) and (abs(ver[i][k].Eta()) < 2.4) and (abs(abs(ver[i][j].Eta()) - abs(ver[i][k].Eta())) < 1.6):
                                   num += 1
                                   break
                         break
     return (num / len(ver))
def doubleJetMassEfficVal(inputlist, ptVal):
     num = 0
     ver = inputlist
     for i in range(len(ver)):
          for j in range(len(ver[i])):
               if (ver[i][j].Pt() >= ptVal) and (abs(ver[i][j].Eta()) < 2.5):
                    for k in range(len(ver[i])):
                         if (ver[i][k].Pt() >= ptVal) and (abs(ver[i][k].Eta()) < 2.5) and (abs(abs(ver[i][j].Eta()) - abs(ver[i][k].Eta())) < 1.5) and ((ver[i][j].M() + ver[i][k].M()) > 300) and (k!=j):
                              num += 1
                              break
                    break
     return (num / len(ver))
def tripleJetEfficVal(inputlist, ptVal):
     num = 0
     ver = inputlist
     for i in range(len(ver)):
          for j in range(len(ver[i])):
               if ver[i][j].Pt() >= 95:
                    for k in range(len(ver[i])):
                         if (ver[i][k].Pt() >= 75) and (abs(ver[i][k].Eta()) < 2.5) and (k!=j):
                              for m in range(len(ver[i])):
                                   if (ver[i][m].Pt() >= ptVal) and (abs(ver[i][k].Eta()) < 2.5) and (m!=k) and (m!=j):
                                        num += 1
                                        break
                              break
                    break
     return (num / len(ver))
def EtmissEfficVal(inputlist, EVal):
     num = 0
     ver = inputlist
     for i in range(len(ver)):
          vectorsumX = np.zeros(1, dtype=object)
          vectorsumY = np.zeros(1, dtype=object)
          for j in range(len(ver[i])):
               if (abs(ver[i][j].Eta()) < 5.0):
                    vectorsumX[0] = vectorsumX[0] + ver[i][j].Px()
                    vectorsumY[0] = vectorsumX[0] + ver[i][j].Px()
          vectorsum = np.sqrt(vectorsumX[0]**2 + vectorsumY[0]**2)
          if (vectorsum > EVal):
               num += 1
     return (num / len(ver))
def HtEfficVal(inputlist, EVal):
     num = 0
     ver = inputlist
     for i in range(len(ver)):
          scalarsum = np.zeros(1, dtype=object)
          for j in range(len(ver[i])):
               if (ver[i][j].Pt() > 30) and (abs(ver[i][j].Eta()) < 2.4):
                    scalarsum[0] = scalarsum[0] + ver[i][j].Pt()
          if (scalarsum[0] > EVal):
               num += 1
     return (num / len(ver))

#Effic Uncertainty Func
def EfficUnc(Npass, Ntotal):
     if (Npass > 0) and (Ntotal > 0):
          unc = np.sqrt((Ntotal**3 + Npass**3) / ((Npass)*(Ntotal**5)))
     else:
          unc = 0
     return unc

def singleJetEfficNpass(inputlist, ptVal): # single jet
     num = 0
     ver = inputlist
     for i in range(len(ver)):
          if len(ver[i]) > 0:
               for j in range(len(ver[i])):
                    if (ver[i][j].Pt() > ptVal) and (abs(ver[i][j].Eta()) < 2.4):
                         num += 1
                         break
     return num
def singleJetEfficNtotal(inputlist, ptVal):
     num = 0
     ver = inputlist
     for i in range(len(ver)):
          if len(ver[i]) > 0:
               for j in range(len(ver[i])):
                    if (ver[i][j].Pt() > ptVal) and (abs(ver[i][j].Eta()) < 2.4):
                         num += 1
                         break
     return len(ver)

def doubleJetEfficNpass(inputlist, ptVal): # double jet
     num = 0
     ver = inputlist
     for i in range(len(ver)):
          if len(ver[i]) >= 2:
               for j in range(len(ver[i])):
                    if (ver[i][j].Pt() > ptVal) and (abs(ver[i][j].Eta()) < 2.4):
                         for k in range(len(ver[i])):
                              if (k!=j) and (ver[i][k].Pt() > ptVal) and (abs(ver[i][k].Eta()) < 2.4) and (abs(abs(ver[i][j].Eta()) - abs(ver[i][k].Eta())) < 1.6):
                                   num += 1
                                   break
                         break
     return num
def doubleJetEfficNtotal(inputlist, ptVal):
     num = 0
     ver = inputlist
     for i in range(len(ver)):
          if len(ver[i]) >= 2:
               for j in range(len(ver[i])):
                    if (ver[i][j].Pt() > ptVal) and (abs(ver[i][j].Eta()) < 2.4):
                         for k in range(len(ver[i])):
                              if (k!=j) and (ver[i][k].Pt() > ptVal) and (abs(ver[i][k].Eta()) < 2.4) and (abs(abs(ver[i][j].Eta()) - abs(ver[i][k].Eta())) < 1.6):
                                   num += 1
                                   break
                         break
     return len(ver)

def doubleJetMassEfficNpass(inputlist, ptVal): # double jet mass
     num = 0
     ver = inputlist
     for i in range(len(ver)):
          for j in range(len(ver[i])):
               if (ver[i][j].Pt() >= ptVal) and (abs(ver[i][j].Eta()) < 2.5):
                    for k in range(len(ver[i])):
                         if (ver[i][k].Pt() >= ptVal) and (abs(ver[i][k].Eta()) < 2.5) and (abs(abs(ver[i][j].Eta()) - abs(ver[i][k].Eta())) < 1.5) and ((ver[i][j].M() + ver[i][k].M()) > 300) and (k!=j):
                              num += 1
                              break
                    break
     return num
def doubleJetMassEfficNtotal(inputlist, ptVal):
     num = 0
     ver = inputlist
     for i in range(len(ver)):
          for j in range(len(ver[i])):
               if (ver[i][j].Pt() >= ptVal) and (abs(ver[i][j].Eta()) < 2.5):
                    for k in range(len(ver[i])):
                         if (ver[i][k].Pt() >= ptVal) and (abs(ver[i][k].Eta()) < 2.5) and (abs(abs(ver[i][j].Eta()) - abs(ver[i][k].Eta())) < 1.5) and ((ver[i][j].M() + ver[i][k].M()) > 300) and (k!=j):
                              num += 1
                              break
                    break
     return len(ver)

def tripleJetEfficNpass(inputlist, ptVal): # triple jet
     num = 0
     ver = inputlist
     for i in range(len(ver)):
          for j in range(len(ver[i])):
               if ver[i][j].Pt() >= 95:
                    for k in range(len(ver[i])):
                         if (ver[i][k].Pt() >= 75) and (abs(ver[i][k].Eta()) < 2.5) and (k!=j):
                              for m in range(len(ver[i])):
                                   if (ver[i][m].Pt() >= ptVal) and (abs(ver[i][k].Eta()) < 2.5) and (m!=k) and (m!=j):
                                        num += 1
                                        break
                              break
                    break
     return num
def tripleJetEfficNtotal(inputlist, ptVal):
     num = 0
     ver = inputlist
     for i in range(len(ver)):
          for j in range(len(ver[i])):
               if ver[i][j].Pt() >= 95:
                    for k in range(len(ver[i])):
                         if (ver[i][k].Pt() >= 75) and (abs(ver[i][k].Eta()) < 2.5) and (k!=j):
                              for m in range(len(ver[i])):
                                   if (ver[i][m].Pt() >= ptVal) and (abs(ver[i][k].Eta()) < 2.5) and (m!=k) and (m!=j):
                                        num += 1
                                        break
                              break
                    break
     return len(ver)

def EtmissEfficNpass(inputlist, EVal): # Etmiss
     num = 0
     ver = inputlist
     for i in range(len(ver)):
          vectorsumX = np.zeros(1, dtype=object)
          vectorsumY = np.zeros(1, dtype=object)
          for j in range(len(ver[i])):
               if (abs(ver[i][j].Eta()) < 5.0):
                    vectorsumX[0] = vectorsumX[0] + ver[i][j].Px()
                    vectorsumY[0] = vectorsumX[0] + ver[i][j].Px()
          vectorsum = np.sqrt(vectorsumX[0]**2 + vectorsumY[0]**2)
          if (vectorsum > EVal):
               num += 1
     return num
def EtmissEfficNtotal(inputlist, EVal):
     num = 0
     ver = inputlist
     for i in range(len(ver)):
          vectorsumX = np.zeros(1, dtype=object)
          vectorsumY = np.zeros(1, dtype=object)
          for j in range(len(ver[i])):
               if (abs(ver[i][j].Eta()) < 5.0):
                    vectorsumX[0] = vectorsumX[0] + ver[i][j].Px()
                    vectorsumY[0] = vectorsumX[0] + ver[i][j].Px()
          vectorsum = np.sqrt(vectorsumX[0]**2 + vectorsumY[0]**2)
          if (vectorsum > EVal):
               num += 1
     return len(ver)

def HtEfficNpass(inputlist, EVal): # Ht
     num = 0
     ver = inputlist
     for i in range(len(ver)):
          scalarsum = np.zeros(1, dtype=object)
          for j in range(len(ver[i])):
               if (ver[i][j].Pt() > 30) and (abs(ver[i][j].Eta()) < 2.4):
                    scalarsum[0] = scalarsum[0] + ver[i][j].Pt()
          if (scalarsum[0] > EVal):
               num += 1
     return num
def HtEfficNtotal(inputlist, EVal):
     num = 0
     ver = inputlist
     for i in range(len(ver)):
          scalarsum = np.zeros(1, dtype=object)
          for j in range(len(ver[i])):
               if (ver[i][j].Pt() > 30) and (abs(ver[i][j].Eta()) < 2.4):
                    scalarsum[0] = scalarsum[0] + ver[i][j].Pt()
          if (scalarsum[0] > EVal):
               num += 1
     return len(ver)

###################################################

def main(args):

    inFileName = args.inFileName
    print("Reading from " + inFileName)

    inFile = r.TFile.Open(inFileName, "READ")

    tree = inFile.Get("ntuple0/objects")
    ver = inFile.Get("ntuple0/objects/vz")

    # Load variables based on inputs and initialize lists to be used later
    eventNum = tree.GetEntries()
    jetNum = 0
    signalPartCount = 0
    jetPartList = []
    jetFullData = []
    trainingFullData = []

    # Scaling Phi for particles relative to their jet
    def signedDeltaPhi(phi1, phi2):
        dPhi = phi1 - phi2
        if dPhi < -np.pi:
            dPhi = 2 * np.pi + dPhi
        elif dPhi > np.pi:
            dPhi = -2 * np.pi + dPhi
        return dPhi

    jetPartsArray = []
    jetDataArray = []
    signalPartArray = []
    missedSignalPartArray = []
    partType = []
    eventjets = []

##########################################
    print("Beginning Jet Construction")
    h_LeadPt =      r.TH1F(name="h_LeadPt", title='h_leadPT', nbinsx=300, xlow=-100, xup=2200) # Pt Plots
    h_SubLeadPt =   r.TH1F(name="h_SubLeadPt", title='h_SubleadPT', nbinsx=300, xlow=-100, xup=2200) # Pt Plots
    h_AllPt =       r.TH1F(name="h_AllJetPt", title='h_AllJetPT', nbinsx=300, xlow=-100, xup=2200) # Pt Plots
    h_LeadPhi =     r.TH1F(name="h_LeadPhi", title='h_LeadPhi', nbinsx=300, xlow=-3.5, xup=3.5) # Phi Plots
    h_SubLeadPhi =  r.TH1F(name="h_SubLeadPhi", title='h_SubLeadPhi', nbinsx=300, xlow=-3.5, xup=3.5) # Phi Plots
    h_AllPhi =      r.TH1F(name="h_AllPhi", title='h_AllPhi', nbinsx=300, xlow=-3.5, xup=3.5) # Phi Plots
    h_LeadEta =     r.TH1F(name="h_LeadEta", title='h_LeadEta', nbinsx=300, xlow=-5, xup=5) # Phi Plots
    h_SubLeadEta =  r.TH1F(name="h_SubLeadEta", title='h_SubLeadEta', nbinsx=300, xlow=-5, xup=5) # Phi Plots
    h_AllEta =      r.TH1F(name="h_AllEta", title='h_AllEta', nbinsx=300, xlow=-5, xup=5) # Phi Plots
    h_LeadMass =    r.TH1F(name="h_LeadMass", title='h_LeadMass', nbinsx=300, xlow=0, xup=500) # Mass Plots
    h_SubLeadMass = r.TH1F(name="h_SubLeadMass", title='h_SubLeadMass', nbinsx=300, xlow=0, xup=500) # Mass Plots
    h_AllMass =     r.TH1F(name="h_AllMass", title='h_AllMass', nbinsx=300, xlow=0, xup=500) # Mass Plots
    start = time.time()

   
    numofevents = args.numofevents
    pbar = tqdm.tqdm(range(int(numofevents)))
#    pbar = tqdm.tqdm(range(int(eventNum)))    
    missedSignalParts = 0
    signalParts = 0
    for entryNum in pbar:
       # pbar.set_description("Jets: " + str(len(jetPartsArray)) + "; Signal Jets: " + str(signalPartCount))
        tree.GetEntry(entryNum)
        ver = tree.vz
        jetlist = []  

        # Loading particle candidates based on PF or PUPPI input
        if not args.usePuppi:
            obj = tree.pf
            verPf = tree.pf_vz
            verPfX = tree.pf_vx
            verPfY = tree.pf_vy
        else:
            obj = tree.pup
            verPf = tree.pup_vz
            verPfX = tree.pup_vx
            verPfY = tree.pup_vy
        jetNum = 0
        bannedParts = []  # List of indices of particles that have already been used by previous jets
        bannedSignalParts = []  # Same deal but with indices within the gen tree corresponding to signal gen particle
  

        # Loops through pf/pup candidates
        for i in range(len(obj)):
            partType.append(obj[i][1]) #adding particle type
            jetPartList = []
            seedParticle = []
            if jetNum >= N_JET_MAX:  # Limited to 12 jets per event at maximum
                jetNum = 0
                break
            if i not in bannedParts:  # Identifies highest avaiable pT particle to use as seed
                tempTLV = obj[i][0]  # Takes TLorentzVector of seed particle to use for jet reconstruction
                if obj[i][1] in [22, 130]:
                    seedParticle.extend(
                        [
                            0.0,
                            verPfX[i],
                            verPfY[i],
                            obj[i][0].Pt(),
                            obj[i][0].Eta(),
                            obj[i][0].Phi(),
                        ]
                    )  # Add in dZ, dX, dY, Particle Pt, Eta, & Phi, last 3 features to be scaled later
                else:
                    seedParticle.extend(
                        [
                            ver[0] - verPf[i],
                            verPfX[i],
                            verPfY[i],
                            obj[i][0].Pt(),
                            obj[i][0].Eta(),
                            obj[i][0].Phi(),
                        ]
                    )  # Add in dZ, dX, dY, Particle Pt, Eta, & Phi, last 3 features to be scaled later
                jetPartList.extend(seedParticle)  # Add particle features to particle list
                bannedParts.append(i)  # Mark this particle as unavailable for other jets
                for j in range(len(obj)):
                    partFts = []
                    if (
                        obj[i][0].DeltaR(obj[j][0]) <= DELTA_R_MATCH and i != j and (j not in bannedParts)
                    ):  # Look for available particles within deltaR<0.4 of seed particle
                        tempTLV = tempTLV + obj[j][0]  # Add to tempTLV
                        if obj[j][1] == 22 or obj[j][1] == 130:
                            partFts.extend(
                                [
                                    0.0,
                                    verPfX[j],
                                    verPfY[j],
                                    obj[j][0].Pt(),
                                    obj[j][0].Eta(),
                                    obj[j][0].Phi(),
                                ]
                            )  # Add in dZ, dX, dY, Particle Pt, Eta, & Phi, last 3 features to be scaled later
                        else:
                            partFts.extend(
                                [
                                    ver[0] - verPf[j],
                                    verPfX[j],
                                    verPfY[j],
                                    obj[j][0].Pt(),
                                    obj[j][0].Eta(),
                                    obj[j][0].Phi(),
                                ]
                            )
                        jetPartList.extend(partFts)  # Add particle features to particle list
                        bannedParts.append(j)  # Mark this particle as unavailable for other jets
               # Ensure all inputs are same length
                while len(jetPartList) < N_PART_PER_JET * N_FEAT:
                    jetPartList.append(0)
                # Store particle inputs and jet features in overall list
                jetPartsArray.append(jetPartList)
                jetDataArray.append((tempTLV.Pt(), tempTLV.Eta(), tempTLV.Phi(), tempTLV.M(), jetPartList[-1]))
           
                h_AllPt.Fill(tempTLV.Pt())
                h_AllPhi.Fill(tempTLV.Phi())
                h_AllEta.Fill(tempTLV.Eta())
                h_AllMass.Fill(tempTLV.M())
                jetlist.append(tempTLV)
                jetNum +=1

        if (len(jetlist)>0):
             h_LeadPt.Fill(jetlist[0].Pt())
             h_LeadPhi.Fill(jetlist[0].Phi())
             h_LeadEta.Fill(jetlist[0].Eta())
             h_LeadMass.Fill(jetlist[0].M())
        if (len(jetlist)>1):
           h_SubLeadPt.Fill(jetlist[1].Pt())
           h_SubLeadPhi.Fill(jetlist[1].Phi())
           h_SubLeadEta.Fill(jetlist[1].Eta())
           h_SubLeadMass.Fill(jetlist[1].M())
        eventjets.append(jetlist)            

    c = r.TCanvas()

    h_LeadPhi.SetLineColor(r.kRed)
    h_LeadPhi.SetTitle("Lead Jet #phi")
    h_LeadPhi.Draw()
    c.Draw()
#    c.SaveAs('h_LeadPhi.png')
    c.Clear()

    h_SubLeadPhi.SetLineColor(r.kGreen)
    h_SubLeadPhi.SetTitle("Sub-Lead Jet #phi")
    h_SubLeadPhi.Draw()
    c.Draw()
#    c.SaveAs('h_SubLeadPhi.png')
    c.Clear()

    h_AllPhi.SetLineColor(r.kBlue)
    h_AllPhi.SetTitle("All Jets #phi")
    h_AllPhi.Draw()
    c.Draw()
#    c.SaveAs('h_AllPhi.png')
    c.Clear()

    h_LeadEta.SetLineColor(r.kRed)
    h_LeadEta.SetTitle("Lead Jet #eta")
    h_LeadEta.Draw()
    c.Draw()
#    c.SaveAs('h_LeadEta.png')
    c.Clear()

    h_SubLeadEta.SetLineColor(r.kGreen)
    h_SubLeadEta.SetTitle("Sub-Lead Jet #eta")
    h_SubLeadEta.Draw()
    c.Draw()
#    c.SaveAs('h_SubLeadEta.png')
    c.Clear()

    h_AllEta.SetLineColor(r.kBlue)
    h_AllEta.SetTitle("All Jets #eta")
    h_AllEta.Draw()
    c.Draw()
#    c.SaveAs('h_AllEta.png')
    c.Clear()

    c.SetLogy()

    h_LeadPt.SetLineColor(r.kRed)
    h_LeadPt.SetTitle("Lead Jet p_{T}")
    h_LeadPt.Draw()
    c.Draw()
#    c.SaveAs('h_LeadPt.png')
    c.Clear()

    h_SubLeadPt.SetLineColor(r.kGreen)
    h_SubLeadPt.SetTitle("Sub-Lead Jet p_{T}")
    h_SubLeadPt.Draw()
    c.Draw()
#    c.SaveAs('h_SubLeadPt.png')
    c.Clear()

    h_AllPt.SetLineColor(r.kBlue)
    h_AllPt.SetTitle("All Jets p_{T}")
    h_AllPt.Draw()
    c.Draw()
#    c.SaveAs('h_AllPt.png')
    c.Clear()
   
    h_LeadMass.SetLineColor(r.kRed)
    h_LeadMass.SetTitle("Lead Jet Mass")
    h_LeadMass.Draw()
    c.Draw()
#    c.SaveAs('h_LeadMass.png')
    c.Clear()

    h_SubLeadMass.SetLineColor(r.kGreen)
    h_SubLeadMass.SetTitle("Sub-Lead Jet Mass")
    h_SubLeadMass.Draw()
    c.Draw()
#    c.SaveAs('h_SubLeadMass.png')
    c.Clear()

    h_AllMass.SetLineColor(r.kBlue)
    h_AllMass.SetTitle("All Jets Mass")
    h_AllMass.Draw()
    c.Draw()
#    c.SaveAs('h_AllMass.png')
    c.Clear()


  
    end = time.time() # timing
    print(f'Elapsed Time: {(end - start)}')

    print(f'Length of eventjets = {len(eventjets)}')

#    print(f'Length of nested list = {len(eventjets[0])}')
#    for i in range(len(eventjets)):
#         if len(eventjets) <= 2:
#              print(f'Event {i} has less than or equal to 2 jets')


    ## Calling Effic Functions
    singleJetEffic(eventjets, 180, 2.4)
    doubleJetEffic(eventjets, 150, 2.5)
    doubleJetdeltaEtaEffic(eventjets, 112, 2.4, 1.6)
    doubleJetMassEffic(eventjets, 160, 35, 620, 5)
    doubleJetMass2Effic(eventjets, 30, 2.5, 1.5, 300)
    tripleJetEffic(eventjets, 95, 75, 65, 2.5)
    EtmissEffic(eventjets, 200, 5.0)
    HtEffic(eventjets, 450, 30, 2.4)
    EtEffic(eventjets, 2000, 5.0)
    
    ## Plotting Effic Curves
    from array import array

    begin = time.time()
    c1 = r.TCanvas('c1', 'Title', 200, 10, 700, 500)
    c1.SetGrid()

    n = 40 # 40
    x1, x2, x3, y1, y2, y3, y4, y5, y6 = array('f'), array('f'), array('f'), array('f'), array('f'), array('f'), array('f'), array('f'), array('f')
    ex, ey1, ey2, ey3, ey4, ey5, ey6  = array('f'), array('f'), array('f'), array('f'), array('f'), array('f'), array('f')
 
    for i in range(n):
         x1.append(10*i) # Input
         x2.append(7.5*i)
         x3.append(133.333*i) 
         y1.append(singleJetEfficVal(eventjets, x1[i])) # Outpit
         y2.append(doubleJetEfficVal(eventjets, x1[i]))
         y3.append(doubleJetMassEfficVal(eventjets, x1[i]))
         y4.append(tripleJetEfficVal(eventjets, x1[i]))
         y5.append(EtmissEfficVal(eventjets, x2[i]))
         y6.append(HtEfficVal(eventjets, x3[i]))
         ex.append(0) # Error Bars
         ey1.append(EfficUnc(singleJetEfficNpass(eventjets, x1[i]), singleJetEfficNtotal(eventjets, x1[i])))
         ey2.append(EfficUnc(doubleJetEfficNpass(eventjets, x1[i]), doubleJetEfficNtotal(eventjets, x1[i])))
         ey3.append(EfficUnc(doubleJetMassEfficNpass(eventjets, x1[i]), doubleJetMassEfficNtotal(eventjets, x1[i])))
         ey4.append(EfficUnc(tripleJetEfficNpass(eventjets, x1[i]), tripleJetEfficNtotal(eventjets, x1[i])))
         ey5.append(EfficUnc(EtmissEfficNpass(eventjets, x2[i]), EtmissEfficNtotal(eventjets, x2[i])))       
         ey6.append(EfficUnc(HtEfficNpass(eventjets, x3[i]), HtEfficNtotal(eventjets, x3[i])))

    g_singleJet = r.TGraphErrors(n, x1, y1, ex, ey1)
    g_singleJet.SetTitle('Single Jet Effic')
    g_singleJet.SetMarkerColor(2)
    g_singleJet.SetMarkerStyle(5)
    g_singleJet.GetXaxis().SetTitle('Pt [GeV]')
    g_singleJet.GetYaxis().SetTitle('Effic')
    g_singleJet.GetYaxis().SetRangeUser(0, 1.0)
    g_singleJet.Draw()
    c1.Update()
#    c1.SaveAs('singleJetEffic.png')
    c1.Clear()

    g_doubleJet = r.TGraphErrors(n, x1, y2, ex, ey2)
    g_doubleJet.SetTitle('Double Jet Effic')
    g_doubleJet.SetMarkerColor(2)
    g_doubleJet.SetMarkerStyle(5)
    g_doubleJet.GetXaxis().SetTitle('Pt [GeV]')
    g_doubleJet.GetYaxis().SetTitle('Effic')
    g_doubleJet.GetYaxis().SetRangeUser(0, 1.0)
    g_doubleJet.Draw()
    c1.Update()
#    c1.SaveAs('doubleJetEffic.png')
    c1.Clear()

    g_doubleJetMass = r.TGraphErrors(n, x1, y3, ex, ey3)
    g_doubleJetMass.SetTitle('Double Jet + Mass Effic')
    g_doubleJetMass.SetMarkerColor(2)
    g_doubleJetMass.SetMarkerStyle(5)
    g_doubleJetMass.GetXaxis().SetTitle('Pt [GeV]')
    g_doubleJetMass.GetYaxis().SetTitle('Effic')
    g_doubleJetMass.GetYaxis().SetRangeUser(0, 1.0)
    g_doubleJetMass.Draw()
    c1.Update()
#    c1.SaveAs('doubleJetMassEffic.png')
    c1.Clear()

    g_tripleJet = r.TGraphErrors(n, x1, y4, ex, ey4)
    g_tripleJet.SetTitle('Triple Jet Effic')
    g_tripleJet.SetMarkerColor(2)
    g_tripleJet.SetMarkerStyle(5)
    g_tripleJet.GetXaxis().SetTitle('Pt [GeV]')
    g_tripleJet.GetYaxis().SetTitle('Effic')
    g_tripleJet.GetYaxis().SetRangeUser(0, 1.0)
    g_tripleJet.Draw()
    c1.Update()
#    c1.SaveAs('tripleJetEffic.png')
    c1.Clear()

    g_Etmiss = r.TGraphErrors(n, x2, y5, ex, ey5)
    g_Etmiss.SetTitle('Etmiss Effic')
    g_Etmiss.SetMarkerColor(2)
    g_Etmiss.SetMarkerStyle(5)
    g_Etmiss.GetXaxis().SetTitle('Pt [GeV]')
    g_Etmiss.GetYaxis().SetTitle('Effic')
    g_Etmiss.GetYaxis().SetRangeUser(0, 1.0)
    g_Etmiss.Draw()
    c1.Update()
#    c1.SaveAs('EtmissEffic.png')
    c1.Clear()

    g_Ht = r.TGraphErrors(n, x3, y6, ex, ey6)
    g_Ht.SetTitle('Ht Effic')
    g_Ht.SetMarkerColor(2)
    g_Ht.SetMarkerStyle(5)
    g_Ht.GetXaxis().SetTitle('Pt [GeV]')
    g_Ht.GetYaxis().SetTitle('Effic')
    g_Ht.GetYaxis().SetRangeUser(0, 1.0)
    g_Ht.Draw()
    c1.Update()
#    c1.SaveAs('HtEffic.png')
    c1.Clear()




    finish = time.time()
    print(f'Effic Curve Plot Time = {finish - begin}')
    print(f'Total RUN Time = {finish - start}')

###########################################################




###############################################


if __name__ == "__main__":
     parser = argparse.ArgumentParser(description="Process arguments")
     parser.add_argument("inFileName", type=str, help="input ROOT file name")
     parser.add_argument("usePuppi", type=bool, help="candidate type (0 for PF, 1 for PUPPI)")
     parser.add_argument('numofevents', type=int, help='number of input events')
     args = parser.parse_args()

     main(args)
