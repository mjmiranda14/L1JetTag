import h5py
import optparse
import numpy as np
import pathlib
import matplotlib.pyplot as plt
import mplhep as hep

stop_samp = 'signalPartsDatasamp.h5'
stop_bkg = 'missedSignalPartsDatasamp.h5'
QCD_bkg = 'jetDatabkg.h5'
dset = [(stop_samp, 'r', 'stop samp'), (stop_bkg, 'g', 'stop bkg'), (QCD_bkg, 'b', 'QCD bkg')]
for idx, feat in enumerate(['pt', 'eta', 'phi', 'm']):
    for i, d in enumerate(dset):
        h5file_path = d[0]
        array_name = 'Data'
        if i == 2:
            array_name = 'Jet ' + array_name
        with h5py.File(h5file_path, 'r') as h5f:
            print(nph5f[array_name])
            feat_array = h5f[array_name][:,idx]
            n, bins = np.histogram(feat_array, bins=100)
            bin_width = bins[1] - bins[0]
            weights = 1/np.sum(n*bin_width)
            n = n*weights
            plt.stairs(n, edges=bins, fill=True, alpha=0.5, label=d[2]) # density=False, histtype='stepfilled', alpha=0.5, facecolor=d[1], label=d[2])
    plt.xlabel(feat)
    plt.xlabel(feat)
    plt.ylabel('normalized count')
    plt.legend(loc='upper right')
    plt.savefig(feat)
    plt.close()

    
output_path = args.parse
stop = 'trainingDatasamp.h5'
QCD = 'trainingDatabkg.h5'
create_file = uproot.recreate(output_path)
with uproot.open(create_file) as file:
    file['particle_level'] = uproot.newtree({"dz": np.float64,
                                             "dx": np.float64,
                                             "dy": np.float64,
                                             "pt": np.float64,
                                             "eta": np.float64,
                                             "phi": np.float64})
    dset = [(stop, 'r', 'stop'), (QCD, 'b', 'QCD')]
    uproot_dict = {}
    for feat in [('dz',8), ('dx', 9), ('dy', 10), ('pt',11), ('eta', 12), ('phi', 13)]:
        plt.figure(figsize=(10, 8))
        for i, d in enumerate(dset):
            h5file_path = d[0]
            array_name = "Training Data"
            #if i == 2:
            #    array_name = 'Jet ' + array_name
            with h5py.File(h5file_path, 'r') as h5f:
                feat_array = h5f[array_name][:,feat[1]::14]
                n, bins = np.histogram(feat_array, bins=100)
                bin_width = bins[1] - bins[0]
                weights = 1/np.sum(n*bin_width)
                n = n*weights
                n_bins = ak.Array([n,bins])
                uproot_dict[feat[0]] = n_bins
                plt.stairs(n, edges=bins, fill=True, alpha=0.5, label=d[2]) # density=False, histtype='stepfilled', alpha=0.5, facecolor=d[1], label=d[2])
        file['particle_level'].extend(uproot_dict)
        plt.xlabel(feat[0])
        plt.ylabel('normalized count')
        plt.legend(loc='upper right')
        plt.savefig(feat[0] + '_particle_level')
        plt.close()
