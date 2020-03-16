# Please cite: Obšteter et al. (2020)


This folder contains:
1.  `AlphaPart_Simulation.R` script to simulate a two-tier pig breeding programs. It uses AlphaSimR to simulate 20 years of burn-in and 20 years of selection in two programmes:

*   Programme 1 = MaleFlow100: The multiplier uses only nucleus males (25)
*   Programme 2 = MaleFlow20: The multiplier uses nucleus (25) and multiplier males (100)
*   It simulate two polygenic traits. Trait 1 has heritability of 0.25 and is measured in the nucleus as well as the multiplier, while trait 2 has heritability of 0.1 and is measured only in the nucleus. The script creates files PedEvalX.csv (X = 1 or 2, marking the programme) holding the pedigree and individuals' information, DataEvalX.csv holding the genetic trend, and AccuraciesX.csv holding the accuracies of the breeding values.


2.  `AlphaPart_Partition.R` script to partition the simulation results. It uses AlphaPart to partition the estimated and true breeding values of the complete breeding programme or separately for the nucleus and the multiplier. It creates files PartitionPN1.csv and PartitionPN2.csv holding the nucleus / multiplier partitions for Programme 1 and Programme 2.

3.  Essentials folder holds the support files (parameter files, scripts) to run the simulation. The folder has to be in the directory, where the simulation is run from. To estimate the breeding values in the simulation, the user should download renumf90 and blupf90 from Ignacy Misztal's group from the University of Georgia (http://nce.ads.uga.edu/html/projects/programs/) and move it into the ./Essentials folder.
