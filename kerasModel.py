import h5py
import numpy
import tensorflow
from keras.models import Sequential
from numpy import expand_dims
from tensorflow.keras.layers import GRU, Conv1D, Dense

# Load in the datasets for training and compiling the sample weights
with h5py.File("trainingDatatest.h5", "r") as hf:
    dataset = hf["Training Data"][:]
with h5py.File("sampleDatatest.h5", "r") as hf:
    sampleData = hf["Sample Data"][:]

# Separate datasets into inputs and outputs, expand the dimensions of the inputs to be used with Conv1D layers
X = dataset[:, 0 : len(dataset[0]) - 1]
y = dataset[:, len(dataset[0]) - 1]
X = expand_dims(X, axis=-1)

# Establish the sample weights
thebins = numpy.linspace(0, 200, 100)
bkgPts = []
sigPts = []
for i in range(len(sampleData)):
    if y[i] == 1:
        sigPts.append(sampleData[i][0])
    if y[i] == 0:
        bkgPts.append(sampleData[i][0])
bkg_counts, binsbkg = numpy.histogram(bkgPts, bins=thebins)
sig_counts, binssig = numpy.histogram(sigPts, bins=thebins)
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
model = Sequential()
model.add(
    Conv1D(
        filters=50,
        kernel_size=14,
        strides=14,
        activation="relu",
        input_shape=(len(dataset[0]) - 1, 1),
    )
)
model.add(Conv1D(filters=50, kernel_size=1, activation="relu"))
model.add(GRU(100, return_sequences=False))
model.add(Dense(50, activation="relu"))
model.add(Dense(10, activation="relu"))
model.add(Dense(1, activation="sigmoid"))
model.compile(loss="binary_crossentropy", optimizer="adam", metrics=["binary_accuracy"])

# Add in the sample weights, 1-to-1 correspondence with training data
# Sample weight of all signal events being equal to 1
# Sample weight of all background events being equal to the sig/bkg ratio at that jet's pT
weights = []
for i in range(len(sampleData)):
    if y[i] == 1:
        weights.append(1)
    if y[i] == 0:
        jetPt = sampleData[i][0]
        tempPt = int(jetPt / 2)
        if tempPt > 98:
            tempPt = 98
        weights.append(a[tempPt])

# Train the network
callback = tensorflow.keras.callbacks.EarlyStopping(monitor="val_loss", verbose=1, patience=5)
model.fit(
    X,
    y,
    epochs=50,
    batch_size=50,
    verbose=2,
    sample_weight=numpy.asarray(weights),
    validation_split=0.20,
    callbacks=[callback],
)

# Save the network
model.save("L1JetTagModel.h5")
