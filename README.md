# Example of using AlphaPart - R implementation of the method for partitioning genetic trends

Obšteter et al. (2020) Genetics Selection and Evolution https://doi.org/10.1186/s12711-021-00600-x

    @article{LaraEtAl2021,
      title = {AlphaPart - R implementation of the method for partitioning genetic trends},
      author = {Obšteter, Jana and Holl, Justin and Hickey, John M. and Gorjanc, Gregor},
      journal = {Genetics Selection and Evolution},
      year = {2021},
      doi = {https://doi.org/10.1186/s12711-021-00600-x}
    }

AlphaPart code repository: https://github.com/AlphaGenes/AlphaPart

This folder contains the code to simulate the analyse the data.

# 1. `Example session`
The folder contains the scripts to simulate a breeding programme using import and partition the genetic trends. The dataset is analyzed in the example session of the paper. The folder contains:

## 1.1 `Simulation_ExampleSession_Public.R`
A script to simulate a breeding programme with three populations in which we select only males for a correlated trait. It uses AlphaSimR to simulate 10 years of burn-in where we select and mate within populations, and 10 years of import where we selection within populations and import males from population 2 and population 3 into population 1. The script next partitions the genetic trend for the correlated (observed) trait and summarizes the results to quantify the contribution of domestic selection and import in population 1.

## 1.1 `SimulatedPedigree.csv`
Data used in the example session in the paper.

# 2. `MultiTier`
The folder contains the scripts to simulate a two-tier breeding programme and partition the genetic trends. The results of this analysis are presented in the result section of the paper. The folder contains:

## 2.1  `AlphaPart_Simulation.R` 
A script to simulate a two-tier pig breeding programs in which we select animals on an index of two traits. It uses AlphaSimR to simulate 20 years of burn-in where we simulate only nucleus and select on phenotype, and 20 years of selection where we simulate two tiers and select on estimated breeding values. The script simulates two programmes:

*   Programme 1 = MaleFlow100: The multiplier uses only nucleus males (25)
*   Programme 2 = MaleFlow20: The multiplier uses nucleus (25) and multiplier males (100)
*   It simulates two polygenic traits. Trait 1 has heritability of 0.25 and is measured in the nucleus as well as the multiplier, while trait 2 has heritability of 0.1 and is measured only in the nucleus. The script creates files PedEvalX.csv (X = 1 or 2, marking the programme) holding the pedigree and individuals' information, DataEvalX.csv holding the genetic trend, and AccuraciesX.csv holding the accuracies of the breeding values.


## 2.2  `AlphaPart_Partition.R` 
A script to partition the simulation results. It uses AlphaPart to partition the estimated and true breeding values of the complete breeding programme or separately for the nucleus and the multiplier. It creates files PartitionPN1.csv and PartitionPN2.csv holding the nucleus / multiplier partitions for Programme 1 and Programme 2.

## 2.3  `Essentials` 
A folder holds the support files (parameter files, scripts) to run the simulation. The folder has to be in the directory, where the simulation is run from. To estimate the breeding values in the simulation, the user should download renumf90 and blupf90 from Ignacy Misztal's group from the University of Georgia (http://nce.ads.uga.edu/html/projects/programs/) and move it into the ./Essentials folder.
