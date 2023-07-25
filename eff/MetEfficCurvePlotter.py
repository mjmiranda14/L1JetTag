import argparse
import time
import numpy as np
import ROOT as r
import tqdm

r.gROOT.SetBatch(1)

def MetEfficVal(EVal):
     num = 0
     inFileName = args.inFileName
     inFile = r.TFile.Open(inFileName, 'READ')
     tree = inFile.Get('ntuple0/objects')
     ver = inFile.Get('ntuple0/objects/vz')
     eventNum = tree.GetEntries()
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
               if vectorsum > EVal:
                    num += 1
     return (num / eventNum)

def EfficUnc(Npass, Ntotal):
     if (Npass > 0) and (Ntotal > 0):
          unc = np.sqrt((Ntotal**3 + Npass**3) / ((Npass)*(Ntotal**5)))
     else:
          unc = 0
     return unc

def MetEfficNpass(EVal):
     num = 0
     inFileName = args.inFileName
     inFile = r.TFile.Open(inFileName, 'READ')
     tree = inFile.Get('ntuple0/objects')
     ver = inFile.Get('ntuple0/objects/vz')
     eventNum = tree.GetEntries()
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
               if vectorsum > EVal:
                    num += 1
     return num


def main(args):

    inFileName = args.inFileName
    print("Reading from " + inFileName)

    inFile = r.TFile.Open(inFileName, "READ")

    tree = inFile.Get("ntuple0/objects")
    ver = inFile.Get("ntuple0/objects/vz")

    # Load variables based on inputs and initialize lists to be used later
    eventNum = tree.GetEntries()


    from array import array

    print("Beginning Efficiency Curve Plotting")
    c1 = r.TCanvas('c1', 'Title', 200, 10, 700, 500)
    c1.SetGrid()

    n = 40 # num of points in TGraph = 40
    x1, y1, ex, ey1 = array('f'), array('f'), array('f'), array('f')

    for i in tqdm.tqdm(range(n)):
         x1.append(25*i) # Input
         y1.append(MetEfficVal(x1[i]))
         ex.append(0) # Error Bars
         ey1.append(EfficUnc(MetEfficNpass(x1[i]), eventNum))

    g_Met = r.TGraphErrors(n, x1, y1, ex, ey1)
    g_Met.SetTitle('MET Effic')
    g_Met.SetMarkerColor(2)
    g_Met.SetMarkerStyle(5)
    g_Met.GetXaxis().SetTitle('Pt [GeV]')
    g_Met.GetYaxis().SetTitle('Effic')
    g_Met.GetYaxis().SetRangeUser(0, 1.0)
    g_Met.Draw()
    c1.Update()
    c1.SaveAs('MetEffic.png')
 #   c1.Clear()

if __name__ == "__main__":
     parser = argparse.ArgumentParser(description="Process arguments")
     parser.add_argument("inFileName", type=str, help="input ROOT file name")
     args = parser.parse_args()

     main(args)



