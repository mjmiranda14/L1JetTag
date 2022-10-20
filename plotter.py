import h5py
import awkward as ak
import uproot
import argparse
import numpy as np
import pathlib
import matplotlib.pyplot as plt
import mplhep as hep

parser = argparse.ArgumentParser()
parser.add_argument('--QCD-part', action='store', type=str, required=False, help='file path for QCD file particle level features')
parser.add_argument('--stop-part', action='store', type=str, required=False, default=100, help='file path for stop file particle level features')
parser.add_argument('--QCD-bkg', action='store', type=str, required=False, help='file path for QCD background')
parser.add_argument('--stop-samp', action='store', type=str, required=False, help='file path for stop file sample')
parser.add_argument('--stop-bkg', action='store', type=str, required=False, help='file path for stop file background')
parser.add_argument('--output', action='store', type=str, required=True, help='output path')
args = parser.parse_args()

output_path = args.output    #"C:\\Users\\eagle\\Research\\commit\\L1METML\\L1JetTag\\"
output_path = output_path + 'graph.root'

create_file = uproot.writing.writable.recreate(output_path)
#with uproot.open(output_path) as file:
mktree_dict = {}
file_list = []
if args.QCD_part == True:   
    QCD_tree = {"QCD_dz": "var * float64",
                "QCD_dx": "var * float64",
                "QCD_dy": "var * float64",
                "QCD_pt": "var * float64",
                "QCD_eta": "var * float64",
                "QCD_phi": "var * float64"}
    mktree_dict = {**mktree_dict, **QCD_tree}
    file_list.append('QCD')
if args. stop_part == True:
    stop_dict = {"stop_dz": "var * float64",
                 "stop_dx": "var * float64",
                 "stop_dy": "var * float64",
                 "stop_pt": "var * float64",
                 "stop_eta": "var * float64",
                 "stop_phi": "var * float64"}
    mktree_dict = {**mktree_dict, **stop_dict}
    file_list.append('stop')

create_file.mktree('particle_level', mktree_dict)
create_file['particle_level'].show()

for kind in file_list:
    for feat in [('dz',8), ('dx', 9), ('dy', 10), ('pt',11), ('eta', 12), ('phi', 13)]:
        n_bins_conc = ak.Array([])
        if kind == "QCD":
            h5file = args.QCD_part #"trainingDatabkg.h5"
            array_name = "Training Data"
        elif kind == "stop":
            h5file = args.stop_part #"trainingDatabkg.h5"
            array_name = "Training Data"
        with h5py.File(h5file, 'r') as h5f:
            feat_array = h5f[array_name][:,feat[1]::14]
            n, bins = np.histogram(feat_array, bins=100)
            bin_width = bins[1] - bins[0]
            weights = 1/np.sum(n*bin_width)
            n = n*weights
        n_bins = ak.Array([n, bins])
        n_bins_conc = ak.concatenate((n_bins_conc, n_bins),axis=0)
        uproot_dict[kind + "_" + feat[0]] = n_bins
print(uproot_dict.keys())
create_file["particle_level"].extend(uproot_dict)        
        

mktree_dict = {}
file_list = []
if args.QCD_bkg == True:
    QCD_bkg_dict = {"bkg_QCD_pt": "var * float64",
                    "bkg_QCD_eta": "var * float64",
                    "bkg_QCD_phi": "var * float64",
                    "bkg_QCD_m": "var * float64"}
    mktree_dict = {**mktree_dict, **QCD_bkg_dict}
    file_list.append("bkg_QCD")
if args.stop_bkg == True:
    stop_bkg_dict = {"bkg_stop_pt": "var * float64",
                     "bkg_stop_eta": "var * float64",
                     "bkg_stop_phi": "var * float64",
                     "bkg_stop_m": "var * float64"}
    mktree_dict = {**mktree_dict, **stop_bkg_dict}
    file_list.append("bkg_stop")
if args.stop_samp == True:
    stop_samp_dict = {"samp_stop_pt": "var * float64",
                      "samp_stop_eta": "var * float64",
                      "samp_stop_phi": "var * float64",
                      "samp_stop_m": "var * float64"}
    mktree_dict = {**mktree_dict, **stop_samp_dict}
    file_list.append("samp_stop")
    
create_file.mktree('jet_level', mktree_dict)
create_file['jet_level'].show()

uproot_dict = {}
for kind in file_list:
    for idx, feat in enumerate(['pt', 'eta', 'phi', 'm']):
        n_bins_conc = ak.Array([])
        if kind == "bkg_QCD":
            h5file = args.QCD_bkg #"missedSignalPartsDatasamp.h5"
            array_name = "Data"
        elif kind == "bkg_stop":
            h5file = args.stop_bkg #"missedSignalPartsDatasamp.h5"
            array_name = "Data"
        elif kind == "samp_stop":
            h5file = args.stop_samp #"signalPartsDatasamp.h5"
            array_name = "Data"
        with h5py.File(h5file, 'r') as h5f:
            feat_array = h5f[array_name][:,idx]
            n, bins = np.histogram(feat_array, bins=100)
            bin_width = bins[1] - bins[0]
            weights = 1/np.sum(n*bin_width)
            n = n*weights
        n_bins = ak.Array([n, bins])
        n_bins_conc = ak.concatenate((n_bins_conc, n_bins),axis=0)
        uproot_dict[kind + "_" + feat] = n_bins
create_file["jet_level"].extend(uproot_dict)   
