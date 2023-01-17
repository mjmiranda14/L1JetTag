import argparse
import time
import h5py
import numpy as np
import ROOT as r
import tqdm

#inFileName = "root://cmseos.fnal.gov///store/user/mequinna/file4russell_2.root"

SIGNAL_PDG_ID = 1000006
MAX_ETA = 2.3
N_JET_MAX = 12
N_FEAT = 14
N_PART_PER_JET = 10
DELTA_R_MATCH = 0.4

r.gROOT.SetBatch(1)


def main(args):
    #tag = "QCD"
    #ptCut = 200
    #trainPercent = 50
    #usePuppi = 0

    inFileName = args.inFileName
    print("Reading from " + inFileName)

    inFile = r.TFile.Open(inFileName, "READ")

    tree = inFile.Get("ntuple0/objects")
    ver = inFile.Get("ntuple0/objects/vz")

    # Load variables based on inputs and initialize lists to be used later
    eventNum = tree.GetEntries()
    ptCut = args.ptCut
    jetNum = 0
    signalPartCount = 0
    jetPartList = []
    trainArray = []
    testArray = []
    jetFullData = []
    trainingFullData = []

    # One-Hot Encoding for Particle Type
    def scalePartType(a, n):
        if n == 11:
            a.extend((1, 0, 0, 0, 0, 0, 0, 0))  # Electron
        elif n == -11:
            a.extend((0, 1, 0, 0, 0, 0, 0, 0))  # Positron
        elif n == 13:
            a.extend((0, 0, 1, 0, 0, 0, 0, 0))  # Muon
        elif n == -13:
            a.extend((0, 0, 0, 1, 0, 0, 0, 0))  # Anti-Muon
        elif n == 22:
            a.extend((0, 0, 0, 0, 1, 0, 0, 0))  # Photon
        elif n == 130:
            a.extend((0, 0, 0, 0, 0, 1, 0, 0))  # Neutral Meson
        elif n == 211:
            a.extend((0, 0, 0, 0, 0, 0, 1, 0))  # Pion
        elif n == -211:
            a.extend((0, 0, 0, 0, 0, 0, 0, 1))  # Anti-Pion
        else:
            a.extend((0, 0, 0, 0, 0, 0, 0, 0))  # Case for unknown particle

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
        

    print("Beginning Jet Construction")
    start = time.time()
    pbar = tqdm.tqdm(range(eventNum))
    
    missedSignalParts = 0
    signalParts = 0
    for entryNum in pbar:
        pbar.set_description("Jets: " + str(len(jetPartsArray)) + "; Signal Jets: " + str(signalPartCount))
        tree.GetEntry(entryNum)
        ver = tree.vz
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
                scalePartType(seedParticle, abs(obj[i][1]))  # One-Hot Encoding Seed Particle Type
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
                        scalePartType(partFts, obj[j][1])  # One-Hot Encoding Particle Type
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
                    if (
                        len(jetPartList) >= N_PART_PER_JET * N_FEAT
                    ):  # If you reach 10 particles in one jet, break and move on
                        break
                if abs(tempTLV.Pt()) < ptCut:  # Neglect to save jet if it falls below pT Cut
                    break
                # Scaling particle pT, Eta, and Phi based on jet pT, Eta, and Phi
                c = N_PART_PER_JET + 1
                while c < len(jetPartList) - 2:
                    jetPartList[c] = jetPartList[c] / tempTLV.Pt()
                    jetPartList[c + 1] = tempTLV.Eta() - jetPartList[c + 1]
                    tempPhi = jetPartList[c + 2]
                    jetPartList[c + 2] = signedDeltaPhi(tempTLV.Phi(), tempPhi)
                    c += N_FEAT
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
                jetNum += 1

        

        for n in range(len(tree.gen)): 
            
            tlv = 0 
            if (
                    (n not in bannedSignalParts)
                    and (abs(tree.gen[n][1]) == SIGNAL_PDG_ID) 
            ):
                missedSignalParts += 1
                tlv = tree.gen[n][0]
                missedSignalPartArray.append((tlv.Pt(), tlv.Eta(), tlv.Phi(), tlv.M(),
                                        tlv.Px(), tlv.Py(), tlv.Pz() ))

            elif (n in bannedSignalParts):
                tlv = tree.gen[n][0]
                signalPartArray.append(( tlv.Pt(), tlv.Eta(), tlv.Phi(), tlv.M(),
                                        tlv.Px(), tlv.Py(), tlv.Pz() ))
                
                
                
    # Break dataset into training/testing data based on train/test split input
    splitIndex = int(float(args.trainPercent) / 100 * len(jetPartsArray))
    trainArray = jetPartsArray[:splitIndex]
    trainingFullData = jetDataArray[:splitIndex]

    testArray = jetPartsArray[splitIndex:]
    jetFullData = jetDataArray[splitIndex:]

    print("Total Jets " + str(len(jetPartsArray)))
    print("Total No. of Matched Jets: " + str(signalPartCount))
    print("No. of Jets in Training Data: " + str(len(trainArray)))
    print("No. of Jets in Testing Data: " + str(len(testArray)))
    print("No. of missed signal particles during reconstruction: " + str(missedSignalParts))
    print("Total No. of signal particles: " + str(signalParts))

    # Final Check
    print("Debug that everything matches up in length")
    assert len(testArray) == len(jetFullData) and len(trainArray) == len(trainingFullData)

    # Save datasets as h5 files
    
    # Testing Data: Particle Inputs for each jet of Shape [...,141]
    with h5py.File("testingData" + str(args.tag) + ".h5", "w") as hf:
        hf.create_dataset("Testing Data", data=testArray)
    # Jet Data: Jet Features (pT, Eta, Phi, Mass) of each testing data jet of shape [...,4]
    with h5py.File("jetData" + str(args.tag) + ".h5", "w") as hf:
        hf.create_dataset("Jet Data", data=jetFullData)
    # Training Data: Particle Inputs for each jet of Shape [...,141]
    with h5py.File("trainingData" + str(args.tag) + ".h5", "w") as hf:
        hf.create_dataset("Training Data", data=trainArray)
    # Sample Data: Jet Features (pT, Eta, Phi, Mass) of each training data jet of shape [...,4]
    with h5py.File("sampleData" + str(args.tag) + ".h5", "w") as hf: 
        hf.create_dataset("Sample Data", data=trainingFullData)
    with h5py.File("signalPartsData" + str(args.tag) +".h5", "w") as hf: 
        hf.create_dataset("Data", data=signalPartArray)
    with h5py.File("missedSignalPartsData"+ str(args.tag) + ".h5", "w") as hf: 
        hf.create_dataset("Data", data=missedSignalPartArray)
    with h5py.File("ParticleTypes" + str(args.tag) + ".h5", "w") as hf: 
        hf.create_dataset("pdgID", data=partType)
        

    end = time.time()
    print(str(end - start))


if __name__ == "__main__":
     parser = argparse.ArgumentParser(description="Process arguments")
     parser.add_argument("inFileName", type=str, help="input ROOT file name")
     parser.add_argument("tag", type=str, help="data/file tag")
     parser.add_argument("ptCut", type=float, help="pT cut applied to individual jets")
     parser.add_argument("trainPercent", type=int, help="fraction (in perecent) of training data (0-100)")
     parser.add_argument("usePuppi", type=bool, help="candidate type (0 for PF, 1 for PUPPI)")

     args = parser.parse_args()

     main(args)
