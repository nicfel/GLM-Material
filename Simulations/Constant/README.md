# Inference of predictor contributions from phylogenetic trees

To recreate the analyses:
1. Run the matlab file _createMasterStepwise.m_, which will create Master input xml's and will output the true values of the coefficients and predictors into files.
2. Run the Master xml files in Beast2 with Master installed.
3. Run the matlab file _createMascotStepwise.m_, which will create Mascot input xml's.
4. Run the matlab file _createDTAConstant.m_, which will create DTA Beast 1 input xml's.
5. Run the Mascot xml file in Beast2 with Mascot installed and DTA xmls using Beast 1.
6. Copy the Mascot log files into the out folder and the dta log files in the dtaout fodler
7. Run _plotScalers.R_ in R to recreate the plots
