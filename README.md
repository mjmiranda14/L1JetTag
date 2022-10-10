# L1BTag
Scripts and code used for the ongoing Level 1 b quark tagging project for CMS

# Reconstructing Jets:

`$ Python3 DataF.py </path/to/file> (using xrootd or another access mode> QCDpt30 30 50 0`

Order is as follows:

path to file (using xrootd: `root://cmsxrootd.fnal.gov///store/...`)

tag = "QCDpt30" or "Stpt30" in this case.

ptCut = 30 (so, >30 GeV)

trainPercent = 50 (50 % training data)

usePuppi = 0 (0 for pf, 1 for PUPPI)

# Training:

Add paths of the training files resulting from DataForge.py.

# ROC Curves:

Inside ROC.py, add paths of the testing data resulting from the DataForge.py.
