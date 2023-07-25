**PURPOSE OF THIS CODE**: 
Given an ntuple dataset, use a seeded clustering algorithm to create at most 12 jets per event, store these in a larger list of nested lists called `eventjets`, and then preform any of the following: jet feature plotter + trigger efficiency + trigger rate + efficiency curves.

**EXPORTING PLOTS**:
The created plots are saved as png files. One option for exporting them is to use *secure copy* to copy the files or a folder containing the files from the LPC to your local system. 

---
**RUNNING THE CODE:**
- The virtual environment described [here](https://github.com/mjmiranda14/L1JetTag/blob/main/README.md) under "*Setting up Conda and cloning this repo*" must be used with this code

- This code uses `argparse` to pass arguments from the command line and allow for easy customization of the code

- Structure of code: `python JetConstructTriggerAnalysis.py "file.root" [0 for PF, 1 for PUPPI] [jet feat plots - 0 for OFF, 1 for ON] [trigger effic - 0 for OFF, 1 for ON] [trigger rate - 0 for OFF, 1 for ON] [effic curves - 0 for OFF, 1 for ON]`

- Command Line Example: `python JetConstructTriggerAnaylsis.py 'file.root' 1 0 1 1 0`
	- this would run the code over the full dataset using PUPPI particles. It would create and save jet feature plots and efficiency curve plots, but it would not print values for the trigger efficiency and rates. 

- **WARNING**: For some reason using the PF particles of the ntuple in jet construction takes significantly longer and renders the code essentially useless for datasets on the scale of tens of thousands of events. Therefore, please use the PUPPI particle option until this is resolved.

- **NOTE:** the option to RUN over a specific number of events is available by un-commenting the following lines: `numofevents = args.numofevents` and `pbar = tqdm.tqdm(range(int(numofevents))) # argparse option to only RUN over a specific number of events` and `parser.add_argument('numofevents', type=int, help='number of input events')` as well as commenting out `pbar = tqdm.tqdm(range(int(eventNum)))`
	- this will also change the way you run the code from the command line:
	- `python JetConstructTriggerAnaylsis.py 'file.root' 1 10000 0 1 1 0`
	- this would run the same options from the above example except over the first 10000 events of the ntuple
---
**UPDATING TRIGGER REQUIREMENT VALUES**:
- With the updated of Phase II L1 Trigger TDRs, it will become necessary to update the trigger requirements of these and any triggers added in the future

- To make this implementation easier, the triggers are defined as function in the first section of the code before the `main()` section under "*Trigger Effic & Rate Functions*"
	- Trigger requirements are passed as arguments of these functions which are called at the end of the code under the "*Calling Trigger Effic and Rate Functions*"

- For example, the `singleJetTrigger(inputlist, ptVal, etaVal)` is defined this way and currently called with `singleJetTrigger(eventjets, 180, 2.4)` where 180 GeV is the pT threshold and 2.4 is the abs(eta) threshold
	- If one wanted to update this to, pt > 250 GeV and eta < 2.2, one would only need to change the values in this line: `singleJetTrigger(eventjets, 250, 2.2)`
---
**JET FEATURE PLOTTER**:
- Assuming the `jet feat plots` is turned ON, [`TH1F`](https://root.cern.ch/root/htmldoc/guides/users-guide/Histograms.html) histogram plots will be made for jet *pT*, *phi*, *eta*, and *mass*.

- For each feature, 3 plots will be made (for a total of 12 plots):
	1. *All Jet*s: takes a datapoint from every constructed jet of every event for the dataset from the *eventjets* list.
	2. *Lead Jet*: takes a datapoint from the first jet (which is also the highest pT jet due to using a seeded clustering algorithm) of every event 
	3. *SubLead Jet*: takes a datapoint from the second jet of every event

- **WARNING:** This clustering algorithm creates jets based on their pT and delta-R values for a given threshold. So, some events may have 0 or even 1 jet in total.
	- Therefore, you should expect a discrepancy in the number of entries (datapoints) for the lead and sublead plots compared to the total events used.

- the 12 plots are saved as png files and are labeled as follows: 'h_AllPt.png', 'h_LeadPt.png', 'h_SubLeadPt.png'

- One may change aspects of the plots (title, axis title and range, bins, etc.) as normal by editing before the plot and canvas are drawn
---
**TRIGGER EFFICIENCY & RATE**:
- Assuming the `trigg effic` and/or the `trigger rate` is turned ON, values for the efficiency and regular (nominal) rate of select triggers are printed out in the terminal.

- The triggers used are defined as a functions before the `main()` section.
	- in each of these "trigger functions" a counting variable is incremented each time a trigger would "fire" for a particular event
	- the value of this counting variable is then used to calculate the efficiency of the trigger
		- Defined as: (# of passed events) / (total # of events)
	- and to calculate the nominal rate of the trigger
		- Defined as: (# of passed events) / (total # of events) * 40 MHz
		- if this rate is less than or equal to 1 MHz the rate is given in kHz

- Triggers used here (total of 9 triggers):
	- Jet Triggers: *single jet*, *double jet*, *triple jet*
		- the double jet trigger has 4 versions for different requirements such as delta-R and mass. 
	- Energy Sum Triggers: *missing transverse energy (MET)*, *Ht*, *Et*

- After jet construction and the `eventjets` list is populated, these trigger functions are called with a set of arguments corresponding to the various trigger requirements in terms of pT, eta, phi, mass, etc.
	- the values of these requirements can be easily modified in the code
	- and whatever values were used will be printed out with the final value in the terminal
	- for example: 
		- `singleJetTrigger(eventjets, 180, 2.4)` sets the trigger to fire if at least one jet in a given event as a pT > 180 GeV and an absolute value of eta < 2.4
		- `singleJetTrigger(eventjets, 250, 2.0)` sets the trigger to fire if at least one jet in a given event as a pT > 250 GeV and an absolute value of eta < 2.0

---
**EFFICIENCY CURVES**:
- Assuming the `effic curve` is turned ON, [`TGraphErrors`](https://root.cern.ch/root/html534/guides/users-guide/Graphs.html) plots will be made for the trigger efficiency

- Of the 9 triggers defined, only 6 have efficiency curves made here
	- these are: *single jet*,  *double jet*, *double jet + mass*, *triple jet*, *MET*, *HT*

- **Calculating these Efficiencies and Uncertainties**
	- to get the datapoints a separate function for each trigger is defined outside the `main()` section that returns only the efficiency value
	- to get the uncertainties a general `EfficUnc` function is defined based on Poisson Error with arguments `(Npass, Ntotal)`
	- separate functions for each trigger are defined to return a value for `Npass` the number of events the trigger fires for, and a value for `Ntotal` the total number of events over which the trigger is run

- **NOTE:** for each efficiency curve only pT is varied; that is, *effic(pT)*
	- and for triggers with multi-jet requirements, only the last jet pT requirement was varied

- the y-axis range of every plot is set to show [0, 1.0] for ease of comparisons

- the 6 plots are saved as png files and are labeled as follows: 'singleJetEffic.png', 'doubleJetEffic.png', etc.

- One may change aspects of the plots (title, axis title and range, bins, etc.) as normal by editing before the plot and canvas are drawn

- **WARNING**: the error bars are working as intended, but these values are too small to be noticeable on the graph
	- one can easily scale these errors if necessary by modifying the error bar array filling by a scale factor of your choice

---
