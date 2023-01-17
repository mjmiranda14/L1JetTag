import h5py
import awkward as ak
import uproot
import argparse
import numpy as np
import pathlib
import matplotlib.pyplot as plt
import mplhep as hep

parser = argparse.ArgumentParser()
parser.add_argument('--hist-file', action='store', type=str, required=False, help='file path for where the histograms are saved')
parser.add_argument('--jet-feat', action='store', type=str, required=False, help='choose which jetfeatures to plot')
parser.add_argument('--part-feat', action='store', type=str, required=False, help='choose which particle features to plot')
parser.add_argument('--jet-file', action='store', type=str, required=False, choices=['bkg_QCD', 'samp_stop', 'samp_bkg'], help='choose if you want to plot jet features from QCD background, stop sample, or stop background')
parser.add_argument('--part-file', action='store', type=str, required=False, choices=['QCD', 'stop'], help='choose if you want to plot particle features from QCD or stop')
args = parser.parse_args()

input_path = args.hist_file #args.hist_file
part_feat = args.part_feat
part_file = args.part_feat
jet_feat = args.jet_feat
jet_file = args.jet_file
color_list = ['r','b','g','y']
with uproot.open(input_path) as upfile:
    upfile['jet_level'].show()
    tree = upfile['particle_level'].arrays()
    for feat in part_feat:
        for color, file in enumerate(part_file):
            branch_name = file + "_" + feat
            n = tree[branch_name][0]
            bins = tree[branch_name][1]
            plt.stairs(n, edges=bins, fill=True, alpha=0.5, color=color_list[color], label=file) # density=False, histtype='stepfilled', alpha=0.5, facecolor=d[1], label=d[2])
        plt.xlabel(feat)
        plt.xlabel(feat)
        plt.ylabel('normalized count')
        plt.legend(loc='upper right')
        plt.show()
        plt.savefig(feat + "_" + 'particle_level')
        plt.close()
        
    tree = upfile['jet_level'].arrays()
    for feat in jet_feat:
        for file in jet_file:
            branch_name = file + "_" + feat
            n = tree[branch_name][0]
            bins = tree[branch_name][1]
            plt.stairs(n, edges=bins, fill=True, alpha=0.5, label=file) # density=False, histtype='stepfilled', alpha=0.5, facecolor=d[1], label=d[2])
        plt.xlabel(feat)
        plt.xlabel(feat)
        plt.ylabel('normalized count')
        plt.legend(loc='upper right')
        plt.show()
        plt.savefig(feat + '_'+'jet_level')
        plt.close()
