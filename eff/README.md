**PURPOSE OF THIS CODE:** 
Given an ntuple dataset of N events, use a seeded clustering algorithm to create at most 12 jets per event, store these in a larger list of nested lists called `eventjets`, and then preform any of the following: jet feature plotter + trigger efficiency + trigger rate + efficiency curves.

**EXPORTING PLOTS**:
The created plots are saved as png files. One option for exporting them is to use *secure copy* to copy the files or a folder containing the files from the LPC to your local system. 

---
**RUNNING THE CODE:**
- The virtual environment shown [here](https://github.com/ucsd-hep-ex/L1JetTag/tree/main) must be used with this code
- This code uses `argparse` to pass arguments from the command line and allow for easy customization of the code
- Structure of code: `python final.py "file.root" [0 for PF, 1 for PUPPI] [desired number of events] [jet feat plots - 0 for OFF, 1 for ON] [trigger effic - 0 for OFF, 1 for ON] [trigger rate - 0 for OFF, 1 for ON] [effic curves - 0 for OFF, 1 for ON]`
- Command Line Example: `python final.py 'file.root' 0 1000 1 0 0 1`
	-  this would run the code over 1000 events using PF particles. It would create and save jet feature plots and efficiency curve plots, but it would not print values for the trigger efficiency and rates.  

---
**JET FEATURE PLOTTER**:
- Assuming the `jet feat plots` is turned ON, plots will be made for jet *pT*, *phi*, *eta*, and *mass*.
- For each feature, 3 plots will be made (for a total of 12 plots):
	1. *All Jet*s: takes a datapoint from every constructed jet of every event for the dataset from the *eventjets* list.
	2. *Lead Jet*: takes a datapoint from the first jet (which is also the highest pT jet due to using a seeded clustering algorithm) of every event 
	3. *SubLead Jet*: takes a datapoint from the second jet of every event
- **WARNING:**
	- This clustering algorithm creates jets based on their pT and delta-R values for a given threshold. So, some events may have 0 or even 1 jet in total.
		- Therefore, you should expect a discrepancy in the number of entries (datapoints) for the lead and sublead plots compared to the total events used.
- the 12 plots are saved as png files and are labeled as follows: 'h_AllPt.png', 'h_LeadPt.png', 'h_SubLeadPt.png'

---
**TRIGGER EFFICIENCY & RATE**:
- 
