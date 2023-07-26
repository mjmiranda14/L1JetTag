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

## Trigger Effic & Rate Functions
def singleJetTrigger(inputlist, ptVal, etaVal):
     num = 0      # counting variable
     ver = inputlist      # easier to type "ver" than "inputlist"
     for i in range(len(ver)):      # looping over every event stored in your 'eventjets' list
          if len(ver[i]) > 0:      # accounting for events with 0 jets
               for j in range(len(ver[i])):      # looping over every jet in a given event [i]
                    if (ver[i][j].Pt() > ptVal) and (abs(ver[i][j].Eta()) < etaVal):
                         num += 1      # incrementing the counting variable
                         break
     if args.triggefficOn:      # prints out trigger efficiency value if this feature is turned on
          print(f'\nThe effic of single jets with pT > {ptVal} and |\u03B7| < {etaVal} is {num / len(ver)} over {len(ver)} events')
     if args.triggrateOn:      # prints out trigger rate value if this feature is turned on
          if (num/len(ver)*40) > 1:
               res = (f'{(num / len(ver)) * 40} MHz')
          else:      # prints out rate in kHz if sufficiently small
               res = (f'{(num / len(ver)) * 40 * 1000} kHz')
          print(f'The nominal trigger rate of single jets with pT > {ptVal} and |\u03B7| < {etaVal} is {res} with {num} events firing over {len(ver)} events\n')

def doubleJetTrigger(inputlist, ptVal, etaVal, deltaEta):
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
     if args.triggefficOn:
          print(f'\nThe effic of double jets + \u0394\u03B7 with pT > {ptVal} and |\u03B7| < {etaVal} and \u0394\u03B7 < {deltaEta} is {num / len(inputlist)} over {len(inputlist)} events')
     if args.triggrateOn:
          if (num/len(ver)*40) > 1:
               res = (f'{(num / len(ver)) * 40} MHz')
          else:
               res = (f'{(num / len(ver)) * 40 * 1000} kHz')
          print(f'The nominal trigger rate of double jets + \u0394\u03B7 with pT > {ptVal} and |\u03B7| < {etaVal} and \u0394\u03B7 < {deltaEta} is {res} with {num} events firing over {len(ver)} events\n')

def quadJetHtTrigger(inputlist, ptVal1, ptVal2, ptVal3, ptVal4, EVal, ptVal, etaVal):
     num = 0
     ver = inputlist
     for i in range(len(ver)):
          for j in range(len(ver[i])):
               if ver[i][j].Pt() >= ptVal1:
                    for k in range(len(ver[i])):
                         if (ver[i][k].Pt() >= ptVal2) and (abs(ver[i][k].Eta()) < etaVal) and (k!=j):
                              for m in range(len(ver[i])):
                                   if (ver[i][m].Pt() >= ptVal3) and (abs(ver[i][k].Eta()) < etaVal) and (m!=k) and (m!=j):
                                        for n in range(len(ver[i])):
                                             if (ver[i][n].Pt() >= ptVal4) and (abs(ver[i][n].Eta()) < etaVal) and (n!=j) and (n!=k) and (n!=m):
                                                  scalarsum = ver[0][0]
                                                  tmpBool = True
                                                  for a in range(len(ver[i])):
                                                       if (ver[i][a].Pt() > ptVal) and (abs(ver[i][a].Eta()) < etaVal):
                                                            if tmpBool:
                                                                 scalarsum = ver[i][a]
                                                                 tmpBool = False
                                                            else:
                                                                 scalarsum = scalarsum + ver[i][a]
                                                  if (scalarsum.Pt() > EVal):
                                                       num += 1
                                                      # break
                                        break
                              break
                    break
     if args.triggefficOn:
          print(f'\nThe effic of quad jets-Ht with pT > {ptVal1}, {ptVal2}, {ptVal3}, {ptVal4} and energy sum > {EVal}; jets pT > {ptVal} and |\u03B7| < {etaVal} is {num / len(inputlist)} over {len(inputlist)} events')         
     if args.triggrateOn:
          if (num/len(ver)*40) > 1:
               res = (f'{(num / len(ver)) * 40} MHz')
          else:
               res = (f'{(num / len(ver)) * 40 * 1000} kHz')
          print(f'The nominal trigger rate of quad jets-Ht with pT > {ptVal1}, {ptVal2}, {ptVal3}, {ptVal4} and energy sum > {EVal}; jets pT > {ptVal} and |\u03B7| < {etaVal} is {res} with {num} events firing over {len(ver)} events\n')

def HtTrigger(inputlist, EVal, ptVal, etaVal):
     num = 0
     ver = inputlist
     for i in range(len(ver)):
          scalarsum = ver[0][0]
          tmpBool = True
          for j in range(len(ver[i])):
               if (ver[i][j].Pt() > ptVal) and (abs(ver[i][j].Eta()) < etaVal):
                    if tmpBool:
                         scalarsum = ver[i][j]
                         tmpBool = False
                    else:
                         scalarsum = scalarsum + ver[i][j]
          if (scalarsum.Pt() > EVal):
               num += 1
     if args.triggefficOn:
          print(f'\nThe effic of Ht energy sum > {EVal} of jets with pT > {ptVal} and |\u03B7| < {etaVal} is {num / len(inputlist)} over {len(inputlist)} events')
     if args.triggrateOn:
          if (num/len(ver)*40) > 1:
                res = (f'{(num / len(ver)) * 40} MHz')
          else:
                res = (f'{(num / len(ver)) * 40 * 1000} kHz')
          print(f'The nominal trigger rate of Ht energy sum > {EVal} of jets with pT > {ptVal} and |\u03B7| < {etaVal} is {res} with {num} events firing over {len(ver)} events\n')

def MetTrigger(EVal):
     num = 0
     inFileName = args.inFileName
     inFile = r.TFile.Open(inFileName, 'READ')
     tree = inFile.Get('ntuple0/objects')
     ver = inFile.Get('ntuple0/objects/vz')
     eventNum = tree.GetEntries()
     met_values = []
     tree.GetEntry(0)
     ver = tree.pup
     for a in range(eventNum):
          tree.GetEntry(a)
          if len(ver) > 1:
               tmpTLVPx = ver[0][0].Px()
               tmpTLVPy = ver[0][0].Py()
               for b in range(1, len(ver)):
                    tmpTLVPx = tmpTLVPx + ver[b][0].Px()
                    tmpTLVPy = tmpTLVPy + ver[b][0].Py()
               vectorsum = np.sqrt(tmpTLVPx**2 + tmpTLVPy**2)
          met_values.append(vectorsum)
     met_values = np.array(met_values)
     mask = met_values > EVal
     num = len(met_values[mask])
     if args.triggefficOn:
          print(f'\nThe effic of MET energy sum > {EVal} of jets is {num / eventNum} over {eventNum} events')
     if args.triggrateOn:
          if (num/len(ver)*40) > 1:
               res = (f'{(num / eventNum) * 40} MHz')
          else:
               res = (f'{(num / eventNum) * 40 * 1000} kHz')
          print(f'The nominal trigger rate of MET energy sum > {EVal} of jets is {res} with {num} events firing over {eventNum} events\n')

## Effic Curve Plotters
def singleJetEfficVal(inputlist, ptVal):
     etaVal = 2.4
     num = 0
     ver = inputlist
     for i in range(len(ver)):
          if len(ver[i]) > 0:
               for j in range(len(ver[i])):
                    if (ver[i][j].Pt() > ptVal) and (abs(ver[i][j].Eta()) < etaVal):
                         num += 1
                         break
     return (num / len(ver))

def doubleJetEfficVal(inputlist, ptVal):
     etaVal, deltaEta = 2.4, 1.6
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
     return (num / len(ver))

def quadJetHtEfficVal(inputlist, EVal):
     ptVal1, ptVal2, ptVal3, ptVal4, ptVal, etaVal = 70, 55, 40, 40, 30, 2.4 
     num = 0
     ver = inputlist
     for i in range(len(ver)):
          for j in range(len(ver[i])):
               if ver[i][j].Pt() >= ptVal1:
                    for k in range(len(ver[i])):
                         if (ver[i][k].Pt() >= ptVal2) and (abs(ver[i][k].Eta()) < etaVal) and (k!=j):
                              for m in range(len(ver[i])):
                                   if (ver[i][m].Pt() >= ptVal3) and (abs(ver[i][k].Eta()) < etaVal) and (m!=k) and (m!=j):
                                        for n in range(len(ver[i])):
                                             if (ver[i][n].Pt() >= ptVal4) and (abs(ver[i][n].Eta()) < etaVal) and (n!=j) and (n!=k) and (n!=m):
                                                  scalarsum = ver[0][0]
                                                  tmpBool = True
                                                  for a in range(len(ver[i])):
                                                       if (ver[i][a].Pt() > ptVal) and (abs(ver[i][a].Eta()) < etaVal):
                                                            if tmpBool:
                                                                 scalarsum = ver[i][a]
                                                                 tmpBool = False
                                                            else:
                                                                 scalarsum = scalarsum + ver[i][a]
                                                  if (scalarsum.Pt() > EVal):
                                                       num += 1
                                                      # break
                                        break
                              break
                    break
     return (num / len(ver))

def HtEfficVal(inputlist, EVal):
     ptVal, etaVal = 30, 2.4
     num = 0
     ver = inputlist
     for i in range(len(ver)):
          scalarsum = ver[0][0]
          tmpBool = True
          for j in range(len(ver[i])):
               if (ver[i][j].Pt() > ptVal) and (abs(ver[i][j].Eta()) < etaVal):
                    if tmpBool:
                         scalarsum = ver[i][j]
                         tmpBool = False
                    else:
                         scalarsum = scalarsum + ver[i][j]
          if (scalarsum.Pt() > EVal):
               num += 1
     return (num / len(ver))

def MetEfficVal(EVal_bins):
     num = 0
     inFileName = args.inFileName
     inFile = r.TFile.Open(inFileName, 'READ')
     tree = inFile.Get('ntuple0/objects')
     ver = inFile.Get('ntuple0/objects/vz')
     eventNum = tree.GetEntries()
     met_values = []
     tree.GetEntry(0)
     ver = tree.pup
     for a in range(eventNum):
          tree.GetEntry(a)
          if len(ver) > 1:
               tmpTLVPx = ver[0][0].Px()
               tmpTLVPy = ver[0][0].Py()
               for b in range(1, len(ver)):
                    tmpTLVPx = tmpTLVPx + ver[b][0].Px()
                    tmpTLVPy = tmpTLVPy + ver[b][0].Py()
               vectorsum = np.sqrt(tmpTLVPx**2 + tmpTLVPy**2)
          met_values.append(vectorsum)
     met_values = np.array(met_values)
     num = np.digitize(met_values, EVal_bins)
     return (num / eventNum)

## Effic Uncertainty Func
def EfficUnc(Npass, Ntotal):
     try:
          unc = np.sqrt((Ntotal**3 + Npass**3) / ((Npass)*(Ntotal**5)))
     except:
          unc = 0
     return unc

def singleJetEfficNpass(inputlist, ptVal):
     etaVal = 2.4
     num = 0
     ver = inputlist
     for i in range(len(ver)):
          if len(ver[i]) > 0:
               for j in range(len(ver[i])):
                    if (ver[i][j].Pt() > ptVal) and (abs(ver[i][j].Eta()) < etaVal):
                         num += 1
                         break
     return num

def doubleJetEfficNpass(inputlist, ptVal):
     etaVal, deltaEta = 2.4, 1.6
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
     return num

def quadJetHtEfficNpass(inputlist, EVal):
     ptVal1, ptVal2, ptVal3, ptVal4, ptVal, etaVal = 70, 55, 40, 40, 30, 2.4
     num = 0
     ver = inputlist
     for i in range(len(ver)):
          for j in range(len(ver[i])):
               if ver[i][j].Pt() >= ptVal1:
                    for k in range(len(ver[i])):
                         if (ver[i][k].Pt() >= ptVal2) and (abs(ver[i][k].Eta()) < etaVal) and (k!=j):
                              for m in range(len(ver[i])):
                                   if (ver[i][m].Pt() >= ptVal3) and (abs(ver[i][k].Eta()) < etaVal) and (m!=k) and (m!=j):
                                        for n in range(len(ver[i])):
                                             if (ver[i][n].Pt() >= ptVal4) and (abs(ver[i][n].Eta()) < etaVal) and (n!=j) and (n!=k) and (n!=m):
                                                  scalarsum = ver[0][0]
                                                  tmpBool = True
                                                  for a in range(len(ver[i])):
                                                       if (ver[i][a].Pt() > ptVal) and (abs(ver[i][a].Eta()) < etaVal):
                                                            if tmpBool:
                                                                 scalarsum = ver[i][a]
                                                                 tmpBool = False
                                                            else:
                                                                 scalarsum = scalarsum + ver[i][a]
                                                  if (scalarsum.Pt() > EVal):
                                                       num += 1
                                                      # break
                                        break
                              break
                    break
     return num

def HtEfficNpass(inputlist, EVal):
     ptVal, etaVal = 30, 2.4
     num = 0
     ver = inputlist
     for i in range(len(ver)):
          scalarsum = ver[0][0]
          tmpBool = True
          for j in range(len(ver[i])):
               if (ver[i][j].Pt() > ptVal) and (abs(ver[i][j].Eta()) < etaVal):
                    if tmpBool:
                         scalarsum = ver[i][j]
                         tmpBool = False
                    else:
                         scalarsum = scalarsum + ver[i][j]
          if (scalarsum.Pt() > EVal):
               num += 1
     return num

def MetEfficNpass(EVal_bins):
     num = 0
     inFileName = args.inFileName
     inFile = r.TFile.Open(inFileName, 'READ')
     tree = inFile.Get('ntuple0/objects')
     ver = inFile.Get('ntuple0/objects/vz')
     eventNum = tree.GetEntries()
     met_values = []
     tree.GetEntry(0)
     ver = tree.pup
     for a in range(eventNum):
          tree.GetEntry(a)
          if len(ver) > 1:
               tmpTLVPx = ver[0][0].Px()
               tmpTLVPy = ver[0][0].Py()
               for b in range(1, len(ver)):
                    tmpTLVPx = tmpTLVPx + ver[b][0].Px()
                    tmpTLVPy = tmpTLVPy + ver[b][0].Py()
               vectorsum = np.sqrt(tmpTLVPx**2 + tmpTLVPy**2)
          met_values.append(vectorsum)
     met_values = np.array(met_values)
     num = np.digitize(met_values, EVal_bins)
     return num


##################################################################################
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
    print("\nBeginning Jet Construction")

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
   
#    numofevents = args.numofevents
#    pbar = tqdm.tqdm(range(int(numofevents))) # argparse option to only RUN over a specific number of events
    pbar = tqdm.tqdm(range(int(eventNum)))    
    missedSignalParts = 0
    signalParts = 0
    for entryNum in pbar:
       # pbar.set_description("Jets: " + str(len(jetPartsArray)) + "; Signal Jets: " + str(signalPartCount))
        tree.GetEntry(entryNum)
        ver = tree.vz
        jetlist = []  

        # Loading particle candidates based on PF or PUPPI input  
        if args.usePuppi:
            obj = tree.pup
            verPf = tree.pup_vz
            verPfX = tree.pup_vx
            verPfY = tree.pup_vy
        else:
            obj = tree.pf
            verPf = tree.pf_vz
            verPfX = tree.pf_vx
            verPfY = tree.pf_vy
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


    end = time.time() # timing for jet construction


#    print(f'Length of nested list = {len(eventjets[0])}')
#    for i in range(len(eventjets)):
#         if len(eventjets) <= 2:
#              print(f'Event {i} has less than or equal to 2 jets')

    begin = time.time()

    ## Plotting Jet Feats
    if args.featplotOn:
         print("Beginning Jet Feature Plotting") 
         c = r.TCanvas()

         h_LeadPhi.SetLineColor(r.kRed)
         h_LeadPhi.SetTitle("Lead Jet #phi")
         h_LeadPhi.Draw()
         c.Draw()
         c.SaveAs('h_LeadPhi.png')
         c.Clear()

         h_SubLeadPhi.SetLineColor(r.kGreen)
         h_SubLeadPhi.SetTitle("Sub-Lead Jet #phi")
         h_SubLeadPhi.Draw()
         c.Draw()
         c.SaveAs('h_SubLeadPhi.png')
         c.Clear()

         h_AllPhi.SetLineColor(r.kBlue)
         h_AllPhi.SetTitle("All Jets #phi")
         h_AllPhi.Draw()
         c.Draw()
         c.SaveAs('h_AllPhi.png')
         c.Clear()

         h_LeadEta.SetLineColor(r.kRed)
         h_LeadEta.SetTitle("Lead Jet #eta")
         h_LeadEta.Draw()
         c.Draw()
         c.SaveAs('h_LeadEta.png')
         c.Clear()

         h_SubLeadEta.SetLineColor(r.kGreen)
         h_SubLeadEta.SetTitle("Sub-Lead Jet #eta")
         h_SubLeadEta.Draw()
         c.Draw()
         c.SaveAs('h_SubLeadEta.png')
         c.Clear()

         h_AllEta.SetLineColor(r.kBlue)
         h_AllEta.SetTitle("All Jets #eta")
         h_AllEta.Draw()
         c.Draw()
         c.SaveAs('h_AllEta.png')
         c.Clear()

         c.SetLogy()

         h_LeadPt.SetLineColor(r.kRed)
         h_LeadPt.SetTitle("Lead Jet p_{T}")
         h_LeadPt.Draw()
         c.Draw()
         c.SaveAs('h_LeadPt.png')
         c.Clear()

         h_SubLeadPt.SetLineColor(r.kGreen)
         h_SubLeadPt.SetTitle("Sub-Lead Jet p_{T}")
         h_SubLeadPt.Draw()
         c.Draw()
         c.SaveAs('h_SubLeadPt.png')
         c.Clear()

         h_AllPt.SetLineColor(r.kBlue)
         h_AllPt.SetTitle("All Jets p_{T}")
         h_AllPt.Draw()
         c.Draw()
         c.SaveAs('h_AllPt.png')
         c.Clear()
   
         h_LeadMass.SetLineColor(r.kRed)
         h_LeadMass.SetTitle("Lead Jet Mass")
         h_LeadMass.Draw()
         c.Draw()
         c.SaveAs('h_LeadMass.png')
         c.Clear()

         h_SubLeadMass.SetLineColor(r.kGreen)
         h_SubLeadMass.SetTitle("Sub-Lead Jet Mass")
         h_SubLeadMass.Draw()
         c.Draw()
         c.SaveAs('h_SubLeadMass.png')
         c.Clear()

         h_AllMass.SetLineColor(r.kBlue)
         h_AllMass.SetTitle("All Jets Mass")
         h_AllMass.Draw()
         c.Draw()
         c.SaveAs('h_AllMass.png')
         c.Clear()
    

    ## Plotting Effic Curves 
    if args.efficcurveOn:
         print("\nBeginning Efficiency Curve Plotting")
         c1 = r.TCanvas('c1', 'Title', 200, 10, 700, 500)
         c1.SetGrid()

         n = 40 # num of points in TGraph

         from array import array
         x1, x3, y1, y2, y3, y5, ex, ey1, ey2, ey3, ey5 = array('f'), array('f'), array('f'),array('f'),array('f'),array('f'),array('f'),array('f'),array('f'),array('f'),array('f')

         for i in range(n):
              x1.append(10*i)
              x3.append(37.5*i)
              y1.append(singleJetEfficVal(eventjets, x1[i]))         
              y2.append(doubleJetEfficVal(eventjets, x1[i]))
              y3.append(quadJetHtEfficVal(eventjets, x3[i]))
              y5.append(HtEfficVal(eventjets, x3[i]))
              ex.append(0)
              ey1.append(EfficUnc(singleJetEfficNpass(eventjets, x1[i]), len(eventjets)))
              ey2.append(EfficUnc(doubleJetEfficNpass(eventjets, x1[i]), len(eventjets)))
              ey3.append(EfficUnc(quadJetHtEfficNpass(eventjets, x3[i]), len(eventjets)))
              ey5.append(EfficUnc(HtEfficNpass(eventjets, x3[i]), len(eventjets)))

 
         g_singleJet = r.TGraphErrors(n, x1, y1, ex, ey1)
         g_singleJet.SetTitle('Single Jet Effic')
         g_singleJet.SetMarkerColor(2)
         g_singleJet.SetMarkerStyle(5)
         g_singleJet.GetXaxis().SetTitle('Pt [GeV]')
         g_singleJet.GetYaxis().SetTitle('Effic')
#         g_singleJet.GetYaxis().SetRangeUser(0, 1.0)
         g_singleJet.Draw()
         c1.Update()
         c1.SaveAs('singleJetEffic.png')
         c1.Clear()

         g_doubleJet = r.TGraphErrors(n, x1, y2, ex, ey2)
         g_doubleJet.SetTitle('Double Jet Effic')
         g_doubleJet.SetMarkerColor(2)
         g_doubleJet.SetMarkerStyle(5)
         g_doubleJet.GetXaxis().SetTitle('Pt [GeV]')
         g_doubleJet.GetYaxis().SetTitle('Effic')
#         g_doubleJet.GetYaxis().SetRangeUser(0, 1.0)
         g_doubleJet.Draw()
         c1.Update()
         c1.SaveAs('doubleJetEffic.png')
         c1.Clear()

         g_quadJetHt = r.TGraphErrors(n, x3, y3, ex, ey3)
         g_quadJetHt.SetTitle('Quad Jet-HT Effic')
         g_quadJetHt.SetMarkerColor(2)
         g_quadJetHt.SetMarkerStyle(5)
         g_quadJetHt.GetXaxis().SetTitle('Pt [GeV]')
         g_quadJetHt.GetYaxis().SetTitle('Effic')
#         g_quadJetHt.GetYaxis().SetRangeUser(0, 1.0)
         g_quadJetHt.Draw()
         c1.Update()
         c1.SaveAs('quadJetHtEffic.png')
         c1.Clear()

         g_Ht = r.TGraphErrors(n, x3, y5, ex, ey5)
         g_Ht.SetTitle('HT Effic')
         g_Ht.SetMarkerColor(2)
         g_Ht.SetMarkerStyle(5)
         g_Ht.GetXaxis().SetTitle('Pt [GeV]')
         g_Ht.GetYaxis().SetTitle('Effic')
#         g_Ht.GetYaxis().SetRangeUser(0, 1.0)
         g_Ht.Draw()
         c1.Update()
         c1.SaveAs('HtEffic.png')
         c1.Clear()

         x2 = np.linspace(1, 1000, 40+1)
         y4 = MetEfficVal(x2)
         ex = np.zeros(41)
         ey4 = EfficUnc(MetEfficNpass(x2), eventNum)

         g_Met = r.TGraphErrors(n, x2, y4, ex, ey4)
         g_Met.SetTitle('MET Effic')
         g_Met.SetMarkerColor(2)
         g_Met.SetMarkerStyle(5)
         g_Met.GetXaxis().SetTitle('Pt [GeV]')
         g_Met.GetYaxis().SetTitle('Effic')
#         g_Met.GetYaxis().SetRangeUser(0, 1.0)
         g_Met.Draw()
         c1.Update()
         c1.SaveAs('MetEffic.png')
         c1.Clear()

  
    ## Calling Trigger Effic and Rate Functions
    if args.triggefficOn or args.triggrateOn:
         if args.triggefficOn:
              if not args.triggrateOn:
                   print('\nBeginning Efficiency Calculations')
              elif args.triggrateOn:
                   print('\nBeginning Efficiency and Rate Calculations')
         if args.triggrateOn:
              if not args.triggefficOn:
                   print('\nBeginning Rate Calculations')
         singleJetTrigger(eventjets, 180, 2.4)
         doubleJetTrigger(eventjets, 112, 2.4, 1.6)
         quadJetHtTrigger(eventjets, 70, 55, 40, 40, 400, 30, 2.4)
         HtTrigger(eventjets, 450, 30, 2.4)
         MetTrigger(200)

    finish = time.time()
    print(f'\nLength of eventjets = {len(eventjets)}')
    print(f'\nJet Construction Time: {(end - start)}')
    if args.featplotOn or args.triggefficOn or args.triggrateOn or args.efficcurveOn:
         print(f'Plotters and Calculators Time = {finish - begin}')
    print(f'Total RUN Time = {finish - start}\n')

###############################################
#    print(f'usePuppi val = {bool(args.usePuppi)}') 
#    print(f'featplotOn val = {bool(args.featplotOn)}')
#    print(f'triggefficOn val = {bool(args.triggefficOn)}')    
#    print(f'triggrateOn val = {bool(args.triggrateOn)}')
#    print(f'efficcurveOn val = {bool(args.efficcurveOn)}')

#    for i in range(len(eventjets)):  # checking for events with less than 2 jets as well as 0 pT lead jets
#         if len(eventjets[i]) < 2:
#              print(f'Event {i} has only {len(eventjets[i])} jets.')
#         try:
#              if eventjets[i][0] == 0.0:
#                   print(f'Event {i} has a 0 pT lead jet')
#         except:
#              continue
###############################################


if __name__ == "__main__":
     parser = argparse.ArgumentParser(description="Process arguments")
     parser.add_argument("inFileName", type=str, help="input ROOT file name")
     parser.add_argument("usePuppi", type=int, help="candidate type (0 for PF, 1 for PUPPI)")
#     parser.add_argument('numofevents', type=int, help='number of input events')
     parser.add_argument('featplotOn', type=int, help='jet feat plotter (0 for OFF, 1 for ON)')
     parser.add_argument('triggefficOn', type=int, help='trigger effic calculator (0 for OFF, 1 for ON)')
     parser.add_argument('triggrateOn', type=int, help='trigger rate calculator (0 for OFF, 1 for ON)')
     parser.add_argument('efficcurveOn', type=int, help='jet feat plotter (0 for OFF, 1 for ON)')
     args = parser.parse_args()

     main(args)
