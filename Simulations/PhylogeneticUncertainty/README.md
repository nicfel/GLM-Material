# Inference of predictor contributions from phylogenetic trees

To recreate the analyses:
1. Run the matlab file _createMasterPhylogeneticUncertainty.m_, which will create Master input xml's and will output the true values of the coefficients and predictors into files.
2. Run the Master xml files in Beast2 with Master installed.
3. Run the matlab code _generateSequences.m_ to simulate sequence alignments using Seq-Gen.
4. Run the matlab file _createMascotPhylogeneticUncertainty_, which will create Mascot input xml's.
5. Run the matlab code _makeTreeInferenceXMLs.m_ to create the xmls which also infer the tree
5. Run the Mascot xml files in Beast2 with Mascot installed.
6. Copy the Mascot log files with fixed trees into the out folder and the ones where the trees are inferred into the treeinfout folder
7. Run _plotScalers.R_ in R to recreate the plots
