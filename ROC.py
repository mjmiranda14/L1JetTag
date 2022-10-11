from numpy import loadtxt
from keras.models import load_model
from sklearn.model_selection import train_test_split
import matplotlib
import h5py
from numpy import expand_dims

matplotlib.use("Agg")
import sys, os, numpy
import tensorflow

# Establish the input and output data for each dataset, expanding dimensions for compatability with Conv1D layer input

# First Dataset
#with h5py.File("/data/t3home000/aidandc/testingDataver3HH.h5", "r") as hf:
 #   dataset = hf["Testing Data"][:]
    
with h5py.File("testingDataStop_1_30.h5", "r") as hf:
    dataset = hf["Testing Data"][:]
with h5py.File("testingDataQCD_30.h5", "r") as hf:
    datasetQCD = hf["Testing Data"][:]
    
dataset = np.concatenate((dataset,datasetQCD)) #Stacking datasets on top of another
np.random.shuffle(dataset) #shuffling rows
   

N_PART_PER_JET = 10
N_FEAT = 14
A = dataset[:, 0 : len(dataset[0]) - 1]
b = dataset[:, len(dataset[0]) - 1]
#A = expand_dims(A, axis=3)
A = A.reshape((A.shape[0], N_PART_PER_JET, N_FEAT))

# Second Dataset
#with h5py.File("/data/t3home000/aidandc/testingDataHHThreeV.h5", "r") as hf:
 #   dataset1 = hf["Testing Data"][:]

#A1 = dataset1[:, 0 : len(dataset1[0]) - 1]
#b1 = dataset1[:, len(dataset1[0]) - 1]
#A1 = expand_dims(A1, axis=3)

# Third Dataset
#with h5py.File("/data/t3home000/aidandc/testingDataHHFiveV.h5", "r") as hf:
 #   dataset2 = hf["Testing Data"][:]

#A2 = dataset2[:, 0 : len(dataset2[0]) - 1]
#b2 = dataset2[:, len(dataset2[0]) - 1]
#A2 = expand_dims(A2, axis=3)

from sklearn.metrics import roc_curve
from sklearn.metrics import auc

import matplotlib.pyplot as plt

# Create plot for ROC
plt.figure(1)
plt.plot([0, 1], [0, 1], "k--")

# Load in respective model for the datasets
#model1 = load_model("modelOne.h5")
#model2 = load_model("modelTwo.h5")
#model3 = load_model("modelThree.h5")
model1 = load_model("L1JetTagModel.h5")

# Creating ROC curves based on model predictions for each dataset
Ab_pred_keras = model1.predict(A).ravel()
fpr_Ab, tpr_Ab, thresholds_Ab = roc_curve(b, Ab_pred_keras)
auc_Ab = auc(fpr_Ab, tpr_Ab)
plt.plot(fpr_Ab, tpr_Ab, label="Pf, AUC={:.3f}".format(auc_Ab))

#Bc_pred_keras = model2.predict(A1).ravel()
#fpr_Bc, tpr_Bc, thresholds_Bc = roc_curve(b1, Bc_pred_keras)
#auc_Bc = auc(fpr_Bc, tpr_Bc)
#plt.plot(fpr_Bc, tpr_Bc, label="dZ+dXY 3 Vertex (area={:.3f})".format(auc_Bc))

#Cd_pred_keras = model3.predict(A2).ravel()
#fpr_Cd, tpr_Cd, thresholds_Cd = roc_curve(b2, Cd_pred_keras)
#auc_Cd = auc(fpr_Cd, tpr_Cd)
#plt.plot(fpr_Cd, tpr_Cd, label="dZ+dXY 5 Vertex (area={:.3f})".format(auc_Cd))

# Establish labels and save image
plt.xlabel("False positive rate")
plt.ylabel("True positive rate")
plt.title("L1 LLP Tag Model ROC Curve", fontsize=16)
plt.legend(loc="best")
plt.savefig("ROCCurves.png")
