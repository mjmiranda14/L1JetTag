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
def singleJetEffic(inputlist, thold):
     num = 0
     ver = inputlist
     for i in range(len(ver)):
          if len(ver[i]) > 0:
               for j in range(len(ver[i])):
                    if ver[i][j].Pt() > thold:
                         num += 1
                         break
     print(f'The effic of single jets with pT > {thold} is {num / len(ver)} over {len(ver)} events')

def doubleJetEffic(inputlist, ptVal, etaVal):
     num = 0
     ver = inputlist
     for i in range(len(ver)):
          if len(ver[i]) >= 2:
               for j in range(len(ver[i])):
                    if (ver[i][j].Pt() > ptVal) and (abs(ver[i][j].Eta()) < etaVal):
                         for k in range(len(ver[i])):
                              if (k!=j) and (ver[i][k].Pt() > ptVal) and (abs(ver[i][k].Eta()) < etaVal):
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

def doubleJetMassEffic(inputlist, ptVal1, ptVal2, massVal):
     num = 0
     ver = inputlist
     for i in range(len(ver)):
          for j in range(len(ver[i])):
               if (ver[i][j].Pt() >= ptVal1):
                    for k in range(len(ver[i])):
                         if (ver[i][k].Pt() >= ptVal2) and ((ver[i][j].M2() + ver[i][k].M2()) > massVal) and (k!=j):
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
                         if (ver[i][k].Pt() >= ptVal) and (abs(ver[i][k].Eta()) < etaVal) and (abs(abs(ver[i][j].Eta()) - abs(ver[i][k].Eta())) < deltaEta) and ((ver[i][j].M2() + ver[i][k].M2()) > massVal) and (k!=j):
                              num += 1
                              break
                    break
     print(f'The effic of double jets + mass with pT > {ptVal} and |\u03B7| < {etaVal} and \u0394\u03B7 < {deltaEta} and M_jj > {massVal} is {num / len(inputlist)} over {len(inputlist)} events')

def tripleJet(inputlist, ptVal1, ptVal2, ptVal3, etaVal):
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


#######################################

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
    h_LeadPt =     r.TH1F(name="h_LeadPt", title='h_leadPT', nbinsx=100, xlow=-100, xup=2000) # Pt Plots
    h_SubLeadPt =  r.TH1F(name="h_SubLeadPt", title='h_SubleadPT', nbinsx=100, xlow=-100, xup=2000) # Pt Plots
    h_AllPt =      r.TH1F(name="h_AllJetPt", title='h_AllJetPT', nbinsx=100, xlow=-100, xup=2000) # Pt Plots
    h_LeadPhi =    r.TH1F(name="h_LeadPhi", title='h_LeadPhi', nbinsx=100, xlow=-3.5, xup=3.5) # Phi Plots
    h_SubLeadPhi = r.TH1F(name="h_SubLeadPhi", title='h_SubLeadPhi', nbinsx=100, xlow=-3.5, xup=3.5) # Phi Plots
    h_AllPhi =     r.TH1F(name="h_AllPhi", title='h_AllPhi', nbinsx=100, xlow=-3.5, xup=3.5) # Phi Plots
    h_LeadEta =    r.TH1F(name="h_LeadEta", title='h_LeadEta', nbinsx=100, xlow=-5, xup=5) # Phi Plots
    h_SubLeadEta = r.TH1F(name="h_SubLeadEta", title='h_SubLeadEta', nbinsx=100, xlow=-5, xup=5) # Phi Plots
    h_AllEta =     r.TH1F(name="h_AllEta", title='h_AllEta', nbinsx=100, xlow=-5, xup=5) # Phi Plots
    start = time.time()

   
    numofevents = args.numofevents
    pbar = tqdm.tqdm(range(int(numofevents)))    
    missedSignalParts = 0
    signalParts = 0
    for entryNum in pbar:
        #pbar.set_description("Jets: " + str(len(jetPartsArray)) + "; Signal Jets: " + str(signalPartCount))
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
                # Add in final value to indicate if particle is matched (1) or unmatched (0)
                # to a gen signal particle by looking for gen signal particles within deltaR<0.4 of jet
                jetPartList.append(0)
                for e in range(len(tree.gen)):
                    if (
                        abs(tree.gen[e][1]) == SIGNAL_PDG_ID
                        and (e not in bannedSignalParts)
                        and abs(tree.gen[e][0].Eta()) < MAX_ETA
                    ):
                        if tree.gen[e][0].DeltaR(tempTLV) <= DELTA_R_MATCH:
                            jetPartList[-1] = 1
                            signalPartCount += 1
                            bannedSignalParts.append(e)
                            break
                # Store particle inputs and jet features in overall list
                jetPartsArray.append(jetPartList)
                jetDataArray.append((tempTLV.Pt(), tempTLV.Eta(), tempTLV.Phi(), tempTLV.M(), jetPartList[-1]))
           
                h_AllPt.Fill(tempTLV.Pt())
                h_AllPhi.Fill(tempTLV.Phi())
                h_AllEta.Fill(tempTLV.Eta())
                jetlist.append(tempTLV)
                jetNum +=1

        if (len(jetlist)>0):
             h_LeadPt.Fill(jetlist[0].Pt())
             h_LeadPhi.Fill(jetlist[0].Phi())
             h_LeadEta.Fill(jetlist[0].Eta())
        #print("eventNum: ",entryNum, "     NJets: ",jetNum, "    len(jetlist): ", len(jetlist))
        #print("Jet 0: ", jetlist[0].Pt(),jetlist[0].Eta(), jetlist[0].Phi())
        #print("Jet 1: ", jetlist[1].Pt(),jetlist[1].Eta(), jetlist[1].Phi())        
        #print("Jet 2: ", jetlist[2].Pt(),jetlist[2].Eta(), jetlist[2].Phi())        
        if (len(jetlist)>1):
           h_SubLeadPt.Fill(jetlist[1].Pt())
           h_SubLeadPhi.Fill(jetlist[1].Phi())
           h_SubLeadEta.Fill(jetlist[1].Eta())
           #print(entryNum)
        eventjets.append(jetlist)            

    c = r.TCanvas()

    h_LeadPhi.SetLineColor(r.kBlue)
    h_LeadPhi.SetTitle("Lead Jet #phi")
    h_LeadPhi.Draw()
    c.Draw()
#    c.SaveAs('h_LeadPhi.png')
    c.Clear()

    h_SubLeadPhi.SetLineColor(r.kBlue)
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

    h_LeadEta.SetLineColor(r.kBlue)
    h_LeadEta.SetTitle("Lead Jet #eta")
    h_LeadEta.Draw()
    c.Draw()
#    c.SaveAs('h_LeadEta.png')
    c.Clear()

    h_SubLeadEta.SetLineColor(r.kBlue)
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

    h_LeadPt.SetLineColor(r.kBlue)
    h_LeadPt.SetTitle("Lead Jet p_{T}")
    h_LeadPt.Draw()
    c.Draw()
#    c.SaveAs('h_LeadPt.png')
    c.Clear()

    h_SubLeadPt.SetLineColor(r.kBlue)
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


  
    end = time.time() # timing
    print(f'Elapsed Time: {str(end - start)}')

    print(f'Length of eventjets = {len(eventjets)}')

#    print(f'Length of nested list = {len(eventjets[0])}')
#    for i in range(len(eventjets)):
#         if len(eventjets) <= 2:
#              print(f'Event {i} has less than or equal to 2 jets')

#    print(f'eventjets[0][0].Phi() = {eventjets[0][0].Phi()}')
#    print(f'eventjets[0][1].Phi() = {eventjets[0][1].Phi()}')
#    print(f'eventjets[{len(eventjets)-1}][0].Phi() = {eventjets[len(eventjets)-1][0].Phi()}')
#    print(f'eventjets[{len(eventjets)-1}][1].Phi() = {eventjets[len(eventjets)-1][1].Phi()}')

#    for i in range(len(eventjets)):
#         if i >= (len(eventjets) - 5):
#              print(eventjets[i][1].Phi())

    ## Calling Effic Functions
    singleJetEffic(eventjets, 180)
    doubleJetEffic(eventjets, 150, 2.5)
    doubleJetdeltaEtaEffic(eventjets, 112, 2.3, 1.6)
    doubleJetMassEffic(eventjets, 110, 35, 620)
    doubleJetMass2Effic(eventjets, 30, 2.5, 1.5, 300)
    tripleJet(eventjets, 95, 75, 65, 2.5)

#    for i in range(len(eventjets)):
#         for j in range(len(eventjets[i])):
#              print(f'M = {eventjets[i][j].M()}, M2 = {eventjets[i][j].M2()}')
#         print('\n')


#    for i in range(len(eventjets)):
#         for j in range(len(eventjets[i])):
#              print(f'{eventjets[i][j].Pt()}, {eventjets[i][j].M())}')
#         print('\n')
###############################################


if __name__ == "__main__":
     parser = argparse.ArgumentParser(description="Process arguments")
     parser.add_argument("inFileName", type=str, help="input ROOT file name")
     parser.add_argument("usePuppi", type=bool, help="candidate type (0 for PF, 1 for PUPPI)")
     parser.add_argument('numofevents', type=int, help='number of input events')
     args = parser.parse_args()

     main(args)
