# L1LLPTag
Scripts and code used for the ongoing Level 1 LLP tagging project for CMS
Setting up Conda

# install conda
<pre>
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh

export PATH="$HOME/miniconda3/bin:$PATH"
conda config --set auto_activate_base false

conda env create -f environment.yml
</pre>

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