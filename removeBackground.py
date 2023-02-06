import h5py
import argparse
import numpy as np

def main(args):
    trainFileName = args.trainFileName
    testFileName = args.testFileName
    print("Reading signal from " + trainFileName)
    print("Reading  from " + testFileName)
    

    TrainF = h5py.File(trainFileName, "r")
    TestF = h5py.File(testFileName, "r")
    
    TrainData = TrainF["Training Data"][:]
    TestData = TestF["Testing Data"][:]

    Dict = { "Signal in Training" : TrainData[:, -1] == 1, "Signal in Testing" : \
        TestData[:, -1] == 1 }


    newTrainData = TrainData[Dict["Signal in Training"]] #all signal jets in training data
    newTestData = TestData[Dict["Signal in Testing"]] #all signal jets in testing data

    with h5py.File("newTestData" + "ST30" + ".h5", "w") as hf:
    hf.create_dataset("Testing Data", data=newTestData)
    
    with h5py.File("newTrainData" + str(args.tag) + ".h5", "w") as hf:
    hf.create_dataset("Training Data", data=newTrainData)

if __name__ == "__main__":
     parser = argparse.ArgumentParser(description="Process arguments")
     parser.add_argument("trainFileName", type=str, help="input .h5 file name")
     parser.add_argument("testFileName", type=str, help="input .h5 file name")
     parser.add_argument("tag", type=str, help="Sample type. E.g, ST30: STop with 30 GeV Pt cut")
     args = parser.parse_args()

     main(args)
 
