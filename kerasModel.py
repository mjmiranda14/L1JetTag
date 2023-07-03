import h5py
import argparse
import numpy as np
import tensorflow
import matplotlib.pyplot as plt
from tensorflow.keras.models import Model
from tensorflow.keras.layers import Conv1D, Dense, Flatten, Input, GlobalAveragePooling1D
from dataForgeScripts.dataForge import N_FEAT, N_PART_PER_JET

def main(args):
    
    signalTrainFile = args.SignalTrainFile
    bkgTrainFile = args.BkgTrainFile
    jetData_TrainFile = args.JetData_TrainFile

    print("Reading signal from " + signalTrainFile)
    print("Reading  from " + bkgTrainFile)
    print("Reading  from " + jetData_TrainFile)

    with h5py.File(inFile1, "r") as hf:
        dataset = hf["Training Data"][:]
    with h5py.File(inFile2, "r") as hf:
        datasetQCD = hf["Training Data"][:]
    with h5py.File(inFile3, "r") as hf:
        sampleData = hf["Sample Data"][:]
    
    dataset = np.concatenate((dataset, datasetQCD))#Put datasets on top of one another
#dataset = np.load("AugTrainingDataPt30.npy")
    np.random.shuffle(dataset) #randomize QCD and Stop samples

    
# Separate datasets into inputs and outputs, expand the dimensions of the inputs to be used with Conv1D layers
    X = dataset[:, 0 : len(dataset[0]) - 1]
    y = dataset[:, len(dataset[0]) - 1]
    X = X.reshape((X.shape[0], N_PART_PER_JET, N_FEAT))

    
# Establish the sample weights

    thebins = np.linspace(0, 200, 100)
    bkgPts = []
    sigPts = []
    for i in range(len(sampleData)):
        if y[i] == 1:
            sigPts.append(sampleData[i][0])
        if y[i] == 0:
            bkgPts.append(sampleData[i][0])
    bkg_counts, binsbkg = np.histogram(bkgPts, bins=thebins)
    sig_counts, binssig = np.histogram(sigPts, bins=thebins)
    a = []
    for i in range(len(bkg_counts)):
        tempSig = float(sig_counts[i])
        tempBkg = float(bkg_counts[i])
        if tempBkg != 0:
            a.append(tempSig / tempBkg)
        if tempBkg == 0:
            a.append(0)
# Normalize the sample weights above a certain pT
    for i in range(42, len(a)):
        a[i] = 0.7

    

# Compile the network
    x = inputs = Input(shape=(N_PART_PER_JET, N_FEAT))
    x = Conv1D(
        filters=50,
        kernel_size=4,
        strides=2,
        activation="relu",
    )(x)
    x = Conv1D(filters=50, kernel_size=4,strides=1, activation="relu")(x)
#x = GlobalAveragePooling1D()(x)
    x = Flatten()(x)
    x = Dense(50, activation="relu")(x)
    x = Dense(10, activation="relu")(x)
    outputs = Dense(1, activation="sigmoid")(x)

    model = Model(inputs=inputs, outputs=outputs)

    initial_learning_rate = 0.001
    lr_schedule = tensorflow.keras.optimizers.schedules.ExponentialDecay(
        initial_learning_rate,
        decay_steps=100000,
        decay_rate=0.96,
        staircase=True)

    model.compile(loss="binary_crossentropy", optimizer=tensorflow.keras.optimizers.Adam(
        learning_rate=lr_schedule), metrics=["binary_accuracy"])

# Add in the sample weights, 1-to-1 correspondence with training data
# Sample weight of all signal events being equal to 1
# Sample weight of all background events being equal to the sig/bkg ratio at that jet's pT
#weights = []
#for i in range(len(sampleData)):
     #   if y[i] == 1:
      #      weights.append(1)
       # if y[i] == 0:
        #    jetPt = sampleData[i][0]
         #   tempPt = int(jetPt / 2)
          #  if tempPt > 98:
           #     tempPt = 98
            #weights.append(a[tempPt])

# Train the network
    callback = tensorflow.keras.callbacks.EarlyStopping(monitor="val_loss", verbose=1, patience=100)
    history=model.fit(
        X,
        y,
        epochs=50,
        batch_size=50,
        verbose=2,
        #sample_weight=np.asarray(weights),
        validation_split=0.20,
        callbacks=[callback],
    )
    plt.figure(figsize=(7,5), dpi=120)
    plt.plot(history.history['loss'], label = 'Train')
    plt.plot(history.history['val_loss'], label = 'Validation')
    plt.title('Model Loss', fontsize=25)
    plt.ylabel('loss')
    plt.xlabel('epoch')
    plt.legend(loc='best')
    plt.savefig("modelLoss.png")

# Save the network
    model.save("L1JetTagModel.h5")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process arguments")

    parser.add_argument("SignalTrainFile", type=str, help="input the signal train file")
    parser.add_argument("BkgTrainFile", type=str, help="input the signal train file")
    parser.add_argument("JetData_TrainFile", type=str, help="input the jet data (named sampleData...h5")

    args = parser.parse_args()
    main(args)
