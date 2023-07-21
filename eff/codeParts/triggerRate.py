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


## Trigger Rate Functions
# Jet Triggers
def singleJetTriggRate(inputlist, ptVal, etaVal):
     num = 0
     ver = inputlist
     for i in range(len(ver)):
          if len(ver[i]) > 0:
               for j in range(len(ver[i])):
                    if (ver[i][j].Pt() > ptVal) and (abs(ver[i][j].Eta()) < etaVal):
                         num += 1
                         break
#     res = (f'{(num / len(ver)) * 40} MHz')
     if (num/len(ver)*40) > 1: 
         res = (f'{(num / len(ver)) * 40} MHz')
     else:
         res = (f'{(num / len(ver)) * 40 * 1000} kHz')
     print(f'The nominal trigger rate of single jets with pT > {ptVal} and |\u03B7| < {etaVal} is {res}')

def doubleJetTriggRate(inputlist, ptVal, etaVal):
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
#     res = (f'{(num / len(ver)) * 40} MHz')
     if (num/len(ver)*40) > 1:
         res = (f'{(num / len(ver)) * 40} MHz')
     else:
         res = (f'{(num / len(ver)) * 40 * 1000} kHz')    
     print(f'The nominal trigger rate of double jets with pT > {ptVal} and |\u03B7| < {etaVal} is {res}')

def doubleJetdeltaEtaTriggRate(inputlist, ptVal, etaVal, deltaEta):
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
#     res = (f'{(num / len(ver)) * 40} MHz')
     if (num/len(ver)*40) > 1:
         res = (f'{(num / len(ver)) * 40} MHz')
     else:
         res = (f'{(num / len(ver)) * 40 * 1000} kHz')
     print(f'The nominal trigger rate of double jets + \u0394\u03B7 with pT > {ptVal} and |\u03B7| < {etaVal} and \u0394\u03B7 < {deltaEta} is {res}')

def doubleJetMassTriggRate(inputlist, ptVal1, ptVal2, massVal, etaVal):
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
#     res = (f'{(num / len(ver)) * 40} MHz')
     if (num/len(ver)*40) > 1:
         res = (f'{(num / len(ver)) * 40} MHz')
     else:
         res = (f'{(num / len(ver)) * 40 * 1000} kHz')
     print(f'The nominal trigger rate of double jets + mass with pT > {ptVal1}, {ptVal2}; two jets pT > {ptVal2} and M_jj > {massVal} is {res}')

def doubleJetMass2TriggRate(inputlist, ptVal, etaVal, deltaEta, massVal):
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
#     res = (f'{(num / len(ver)) * 40} MHz')
     if (num/len(ver)*40) > 1:
         res = (f'{(num / len(ver)) * 40} MHz')
     else:
         res = (f'{(num / len(ver)) * 40 * 1000} kHz')
     print(f'The nominal trigger rate of double jets + mass with pT > {ptVal} and |\u03B7| < {etaVal} and \u0394\u03B7 < {deltaEta} and M_jj > {massVal} is {res}')

def tripleJetTriggRate(inputlist, ptVal1, ptVal2, ptVal3, etaVal):
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
#     res = (f'{(num / len(ver)) * 40} MHz')
     if (num/len(ver)*40) > 1:
         res = (f'{(num / len(ver)) * 40} MHz')
     else:
         res = (f'{(num / len(ver)) * 40 * 1000} kHz')
     print(f'The nominal trigger rate of triple jets with pT > {ptVal1}, {ptVal2}, {ptVal3}; two jets pT > {ptVal2}, {ptVal3} and |\u03B7| < {etaVal} is {res}')         

# Energy Sum Triggers
def EtmissTriggRate(inputlist, EVal, etaVal):
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
#     res = (f'{(num / len(ver)) * 40} MHz')
     if (num/len(ver)*40) > 1:
         res = (f'{(num / len(ver)) * 40} MHz')
     else:
         res = (f'{(num / len(ver)) * 40 * 1000} kHz')
     print(f'The nominal trigger rate of E_miss_t energy sum > {EVal} of jets with |\u03B7| < {etaVal} is {res}')

def HtTriggRate(inputlist, EVal, ptVal, etaVal):
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
#     res = (f'{(num / len(ver)) * 40} MHz')
     if (num/len(ver)*40) > 1:
         res = (f'{(num / len(ver)) * 40} MHz')
     else:
         res = (f'{(num / len(ver)) * 40 * 1000} kHz')
     print(f'The nominal trigger rate of H_t energy sum > {EVal} of jets with pT > {ptVal} and |\u03B7| < {etaVal} is {res}')

def EtTriggRate(inputlist, EVal, etaVal):
     num = 0
     ver = inputlist
     for i in range(len(ver)):
          scalarsum = np.zeros(1, dtype=object)
          for j in range(len(ver[i])):
               if (abs(ver[i][j].Eta()) < etaVal):
                    scalarsum[0] = scalarsum[0] + ver[i][j].Pt()
          if (scalarsum[0] > EVal):
               num += 1
#     res = (f'{(num / len(ver)) * 40} MHz')
     if (num/len(ver)*40) > 1:
         res = (f'{(num / len(ver)) * 40} MHz')
     else:
         res = (f'{(num / len(ver)) * 40 * 1000} kHz')
     print(f'The effic of E_t energy sum > {EVal} of jets with |\u03B7| < {etaVal} is {res}')

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
#    h_LeadPt =      r.TH1F(name="h_LeadPt", title='h_leadPT', nbinsx=300, xlow=-100, xup=2200) # Pt Plots
#    h_SubLeadPt =   r.TH1F(name="h_SubLeadPt", title='h_SubleadPT', nbinsx=300, xlow=-100, xup=2200) # Pt Plots
#    h_AllPt =       r.TH1F(name="h_AllJetPt", title='h_AllJetPT', nbinsx=300, xlow=-100, xup=2200) # Pt Plots
#    h_LeadPhi =     r.TH1F(name="h_LeadPhi", title='h_LeadPhi', nbinsx=300, xlow=-3.5, xup=3.5) # Phi Plots
#    h_SubLeadPhi =  r.TH1F(name="h_SubLeadPhi", title='h_SubLeadPhi', nbinsx=300, xlow=-3.5, xup=3.5) # Phi Plots
#    h_AllPhi =      r.TH1F(name="h_AllPhi", title='h_AllPhi', nbinsx=300, xlow=-3.5, xup=3.5) # Phi Plots
#    h_LeadEta =     r.TH1F(name="h_LeadEta", title='h_LeadEta', nbinsx=300, xlow=-5, xup=5) # Phi Plots
#    h_SubLeadEta =  r.TH1F(name="h_SubLeadEta", title='h_SubLeadEta', nbinsx=300, xlow=-5, xup=5) # Phi Plots
#    h_AllEta =      r.TH1F(name="h_AllEta", title='h_AllEta', nbinsx=300, xlow=-5, xup=5) # Phi Plots
#    h_LeadMass =    r.TH1F(name="h_LeadMass", title='h_LeadMass', nbinsx=300, xlow=0, xup=500) # Mass Plots
#    h_SubLeadMass = r.TH1F(name="h_SubLeadMass", title='h_SubLeadMass', nbinsx=300, xlow=0, xup=500) # Mass Plots
#    h_AllMass =     r.TH1F(name="h_AllMass", title='h_AllMass', nbinsx=300, xlow=0, xup=500) # Mass Plots
 

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
           
#                h_AllPt.Fill(tempTLV.Pt())
#                h_AllPhi.Fill(tempTLV.Phi())
#                h_AllEta.Fill(tempTLV.Eta())
#                h_AllMass.Fill(tempTLV.M())
                jetlist.append(tempTLV)
                jetNum +=1

#        if (len(jetlist)>0):
#             h_LeadPt.Fill(jetlist[0].Pt())
#             h_LeadPhi.Fill(jetlist[0].Phi())
#             h_LeadEta.Fill(jetlist[0].Eta())
#             h_LeadMass.Fill(jetlist[0].M())
#        if (len(jetlist)>1):
#           h_SubLeadPt.Fill(jetlist[1].Pt())
#           h_SubLeadPhi.Fill(jetlist[1].Phi())
#           h_SubLeadEta.Fill(jetlist[1].Eta())
#           h_SubLeadMass.Fill(jetlist[1].M())
        eventjets.append(jetlist)            

#    c = r.TCanvas()

#    h_LeadPhi.SetLineColor(r.kRed)
#    h_LeadPhi.SetTitle("Lead Jet #phi")
#    h_LeadPhi.Draw()
#    c.Draw()
#    c.SaveAs('h_LeadPhi.png')
#    c.Clear()

#    h_SubLeadPhi.SetLineColor(r.kGreen)
#    h_SubLeadPhi.SetTitle("Sub-Lead Jet #phi")
#    h_SubLeadPhi.Draw()
#    c.Draw()
#    c.SaveAs('h_SubLeadPhi.png')
#    c.Clear()

#    h_AllPhi.SetLineColor(r.kBlue)
#    h_AllPhi.SetTitle("All Jets #phi")
#    h_AllPhi.Draw()
#    c.Draw()
#    c.SaveAs('h_AllPhi.png')
#    c.Clear()

#    h_LeadEta.SetLineColor(r.kRed)
#    h_LeadEta.SetTitle("Lead Jet #eta")
#    h_LeadEta.Draw()
#    c.Draw()
#    c.SaveAs('h_LeadEta.png')
#    c.Clear()

#    h_SubLeadEta.SetLineColor(r.kGreen)
#    h_SubLeadEta.SetTitle("Sub-Lead Jet #eta")
#    h_SubLeadEta.Draw()
#    c.Draw()
#    c.SaveAs('h_SubLeadEta.png')
#    c.Clear()

#    h_AllEta.SetLineColor(r.kBlue)
#    h_AllEta.SetTitle("All Jets #eta")
#    h_AllEta.Draw()
#    c.Draw()
#    c.SaveAs('h_AllEta.png')
#    c.Clear()

#    c.SetLogy()

#    h_LeadPt.SetLineColor(r.kRed)
#    h_LeadPt.SetTitle("Lead Jet p_{T}")
#    h_LeadPt.Draw()
#    c.Draw()
#    c.SaveAs('h_LeadPt.png')
#    c.Clear()

#    h_SubLeadPt.SetLineColor(r.kGreen)
#    h_SubLeadPt.SetTitle("Sub-Lead Jet p_{T}")
#    h_SubLeadPt.Draw()
#    c.Draw()
#    c.SaveAs('h_SubLeadPt.png')
#    c.Clear()

#    h_AllPt.SetLineColor(r.kBlue)
#    h_AllPt.SetTitle("All Jets p_{T}")
#    h_AllPt.Draw()
#    c.Draw()
#    c.SaveAs('h_AllPt.png')
#    c.Clear()
   
#    h_LeadMass.SetLineColor(r.kRed)
#    h_LeadMass.SetTitle("Lead Jet Mass")
#    h_LeadMass.Draw()
#    c.Draw()
#    c.SaveAs('h_LeadMass.png')
#    c.Clear()

#    h_SubLeadMass.SetLineColor(r.kGreen)
#    h_SubLeadMass.SetTitle("Sub-Lead Jet Mass")
#    h_SubLeadMass.Draw()
#    c.Draw()
#    c.SaveAs('h_SubLeadMass.png')
#    c.Clear()

#    h_AllMass.SetLineColor(r.kBlue)
#    h_AllMass.SetTitle("All Jets Mass")
#    h_AllMass.Draw()
#    c.Draw()
#    c.SaveAs('h_AllMass.png')
#    c.Clear()


  
    end = time.time() # timing
    print(f'Elapsed Time: {(end - start)}')

    print(f'Length of eventjets = {len(eventjets)}')

#    print(f'Length of nested list = {len(eventjets[0])}')
#    for i in range(len(eventjets)):
#         if len(eventjets) <= 2:
#              print(f'Event {i} has less than or equal to 2 jets')


    ## Calling Trigger Rate Functions
    singleJetTriggRate(eventjets, 180, 2.4)
    doubleJetTriggRate(eventjets, 150, 2.5)
    doubleJetdeltaEtaTriggRate(eventjets, 112, 2.4, 1.6)
    doubleJetMassTriggRate(eventjets, 160, 35, 620, 5)
    doubleJetMass2TriggRate(eventjets, 30, 2.5, 1.5, 300)
    tripleJetTriggRate(eventjets, 95, 75, 65, 2.5)
    EtmissTriggRate(eventjets, 200, 5.0)
    HtTriggRate(eventjets, 450, 30, 2.4)
    EtTriggRate(eventjets, 2000, 5.0)
    
###############################################


if __name__ == "__main__":
     parser = argparse.ArgumentParser(description="Process arguments")
     parser.add_argument("inFileName", type=str, help="input ROOT file name")
     parser.add_argument("usePuppi", type=bool, help="candidate type (0 for PF, 1 for PUPPI)")
     parser.add_argument('numofevents', type=int, help='number of input events')
     args = parser.parse_args()

     main(args)
