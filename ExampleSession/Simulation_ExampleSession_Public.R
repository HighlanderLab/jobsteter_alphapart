# 
# ---- Environment ----
# Set working environment
#setwd()

# Load the libraries
library(package = "AlphaSimR")
library(package = "tidyverse")
library(package = "AlphaPart")
library(package = "tidyr")
library(package = "gridExtra")
library(package = "grid")


# ---- General parameters ----
# Set the number of males and females selected in each population
# This simulation includes only selection on males
nMales_pop1   =  5
nFemales_pop1 = 500
nMales_pop2   =  10
nFemales_pop2 = 1000
nMales_pop3   =  10
nFemales_pop3 = 2000
# Set the number of burn-in and evaluation generations
nGenerationBurn = 10
nGenerationEval = 20

# Set the traits to record - we are simulating three correlated traits
GenMeanCols = c("GenMeanT1","GenMeanT2", "GenMeanT3")
GenVarCols  = c("GenVarT1", "GenVarT2", "GenVarT3")

# Simulate the founder population
founderPop = runMacs(nInd = (nMales_pop1+nMales_pop2+nMales_pop3+nFemales_pop1+nFemales_pop2+nFemales_pop3),
                      nChr = 10,
                      segSites = 1000,
                      nThreads = 4,
                      species = "GENERIC")

# Save the founder population
save.image(file = "~/Documents/PhD/Projects/inProgress/AlphaPart/Revision3/SimulationExample/FounderPop_generic.R")

##############################################################################
##############################################################################
# ---- Simulation ----
##############################################################################
##############################################################################
# ---- Simulation/Base population parameters ----
SP = SimParam$new(founderPop)
# Set the genetic variance with a high correlation between traits
VarA = matrix(data = c(1.0, 0.9, 0.8,
                       0.9, 1.0, 0.8,
                       0.8, 0.8, 1.0), nrow = 3); cov2cor(VarA)
# Set the environmental variances to achieve 0.7 and 0.9 h2 for the traits
VarE = matrix(data = c(1.0/0.7-1, 0.0, 0.0, 
                       0.0, 1.0/0.9-1, 0.0,
                       0.0, 0.0, 1.0/0.8-1), nrow = 3); cov2cor(VarE)
# Set the phenotypic variance and check the h2
VarP = VarA + VarE; diag(VarA) / diag(VarP)
# Add an additive trait to the population with the var-cov structure
SP$addTraitA(nQtlPerChr = 1000, mean = c(0, 0, 0), var = diag(VarA), cor = cov2cor(VarA))
# Set gender as systematic (half F - half M)
SP$setGender(gender = "yes_sys")

# ---- Base GN population ----
Base = newPop(founderPop)[1:(nMales_pop1+nMales_pop2+nMales_pop3+nFemales_pop1+nFemales_pop2+nFemales_pop3)]
Base@gender <- sample(c(rep("F", nFemales_pop1+nFemales_pop2+nFemales_pop3), rep("M", nMales_pop1+nMales_pop2+nMales_pop3)))
# Here select the base females and males for both populations - pop1 and pop2
BaseMales_pop1   = Base[Base@gender == "M"][1:nMales_pop1]
BaseFemales_pop1 = Base[Base@gender == "F"][1:nFemales_pop1]
BaseMales_pop2   = Base[Base@gender == "M"][(nMales_pop1+1):(nMales_pop1+nMales_pop2)]
BaseFemales_pop2 = Base[Base@gender == "F"][(nFemales_pop1+1):(nFemales_pop1 + nFemales_pop2)]
BaseMales_pop3   = Base[Base@gender == "M"][(nMales_pop1+nMales_pop2+1):(nMales_pop1 + nMales_pop2 + nMales_pop3)]
BaseFemales_pop3 = Base[Base@gender == "F"][(nFemales_pop1+nFemales_pop2+1):(nFemales_pop1 + nFemales_pop2 + nFemales_pop3)]
# Check the number of individuals in each population
BaseFemales_pop1@nInd
BaseFemales_pop2@nInd
BaseMales_pop1@nInd
BaseMales_pop2@nInd
sum(BaseFemales_pop1@id %in% BaseFemales_pop2@id)
sum(BaseMales_pop1@id %in% BaseMales_pop2@id)
sum(BaseMales_pop2@id %in% BaseMales_pop3@id)
rm(Base)
  
##############################################################################
# ---- Burn-in ----
##############################################################################
# In the burn-in you create three populations - pop1, pop2 and pop3
# Select only the males - pop1 on trait 1, pop2 on trait 2and pop3 on trait 3
DataBurn = data.frame(Generation = rep(1:nGenerationBurn, each=6), 
                      Gender = rep(c("F", "M"), nGenerationBurn*3),
                      Population = rep(c("Pop1", "Pop1", "Pop2", "Pop2", "Pop3", "Pop3"), nGenerationBurn),
                      GenMeanT1 = NA, GenMeanT2 = NA, GenMeanT3 = NA, 
                      GenVarT1  = NA, GenVarT2  = NA, GenVarT3  = NA)

# Initiate a tibble to hold the burn-in information
PedBurnIn <- tibble()

# For each generation of the burn-in perform within population selection and mating
for (Generation in 1:nGenerationBurn) {
  for (pop in c("Pop1", "Pop2", "Pop3")) {
  
    # Set breeding individuals according to the population
    if (pop == "Pop1") {
      BaseFemales = BaseFemales_pop1
      BaseMales = BaseMales_pop1
    } else if (pop == "Pop2") {
      BaseFemales = BaseFemales_pop2
      BaseMales = BaseMales_pop2
    } else if (pop == "Pop3") {
      BaseFemales = BaseFemales_pop3
      BaseMales = BaseMales_pop3
    }
    
    # Mate the individuals
    SelCand = randCross2(females = BaseFemales, males = BaseMales,
                         nCrosses = BaseFemales@nInd, nProgeny = 2)
    
    # Save metrics
    for (gender in c("F", "M")) {
      DataBurn[(DataBurn$Generation == Generation) & (DataBurn$Gender == gender) & (DataBurn$Pop == pop), GenMeanCols] =
        colMeans(SelCand@gv[SelCand@gender == gender,])
      DataBurn[(DataBurn$Generation == Generation) & (DataBurn$Gender == gender) & (DataBurn$Pop == pop), GenVarCols] =
        diag(var(SelCand@gv[SelCand@gender == gender,]))
    }
    
    
    # Phenotype
    SelCand = setPheno(pop = SelCand, varE = VarE)

    # Track the pedigree and related info
    PedBurnIn = rbind(PedBurnIn,
                    tibble(Generation = Generation,
                           IId        = SelCand@id,
                           FId        = SelCand@father,
                           MId        = SelCand@mother,
                           Gender     = SelCand@gender,
                           Program    = "BurnIn",
                           Population = pop,
                           PhenoT1    = SelCand@pheno[,1],
                           PhenoT2    = SelCand@pheno[,2],
                           PhenoT3    = SelCand@pheno[,3],
                           TbvT1      = SelCand@gv[, 1],
                           TbvT2      = SelCand@gv[, 2],
                           TbvT3      = SelCand@gv[, 3]))
    
    # Select new parents according to the population
    if (pop == "Pop1") {
      BaseMales_pop1   = selectInd(pop = SelCand, nInd = nMales_pop1,   gender = "M",
                              use = "pheno", trait = 1)
      BaseFemales_pop1 = selectInd(pop = SelCand, nInd = nFemales_pop1, gender = "F")
    } else if (pop == "Pop2") {
      BaseMales_pop2   = selectInd(pop = SelCand, nInd = nMales_pop2,   gender = "M",
                                   use = "pheno", trait = 2)
      BaseFemales_pop2 = selectInd(pop = SelCand, nInd = nFemales_pop2, gender = "F")
    } else if (pop == "Pop3") {
      BaseMales_pop3   = selectInd(pop = SelCand, nInd = nMales_pop3,   gender = "M",
                                   use = "pheno", trait = 3)
      BaseFemales_pop3 = selectInd(pop = SelCand, nInd = nFemales_pop3, gender = "F")
    }

  }
}

# Plot genetic means
DataBurn %>%
gather(key = "Metric", value = "Value", GenMeanCols) %>%
ggplot(., aes(Generation, Value, color = Metric)) +
  geom_line() + facet_grid(rows = vars(Population), cols = vars(Gender)) + 
  ylab(label = "Genetic mean")

# Plot only the one correlated (observed) trait, i.e. trait 1 for pop2, trait 2 for pop2 and trait 3 for pop3
OneTrait <- DataBurn %>% gather(key = "Metric", value = "Value", GenMeanCols)
rbind(OneTrait[OneTrait$Population == "Pop1" & OneTrait$Metric == "GenMeanT1",], 
      OneTrait[OneTrait$Population == "Pop2" & OneTrait$Metric == "GenMeanT2",],
      OneTrait[OneTrait$Population == "Pop3" & OneTrait$Metric == "GenMeanT3",]) %>% 
ggplot(., aes(Generation, Value, color = Metric)) +
geom_line() + facet_grid(cols = vars(Gender)) + 
ylab(label = "Genetic mean")


# Plot genetic variances
DataBurn %>% gather(key = "Metric", value= "Value", GenVarCols) %>%
ggplot(., aes(Generation, Value, color = Metric)) +
geom_line() +  facet_grid(rows = vars(Population), cols = vars(Gender)) + 
ylab(label = "Genetic variance")


##############################################################################
# ---- Import ----
##############################################################################
# Change the selection scheme to import the genetic material from pop2 and pop3 into pop1

# Initiate the tibble to hold the pedigree
PedImport <- PedBurnIn
# Initiate a dataframe to hold the genetic mean of the import generations
DataImport = data.frame(Generation = rep((nGenerationBurn+1):(nGenerationBurn + nGenerationEval), each=6), 
                      Gender = rep(c("F", "M"), nGenerationEval*3),
                      Population = rep(c("Pop1", "Pop1", "Pop2", "Pop2", "Pop3", "Pop3"), nGenerationBurn),
                      GenMeanT1 = NA, GenMeanT2 = NA, GenMeanT3 = NA, 
                      GenVarT1  = NA, GenVarT2  = NA, GenVarT3  = NA)

# Initiate a dataframe to hold the accuracies
accuracies <- data.frame(Population = NA, Generation = NA, Trait = NA, Cor = NA)

# Set the total import percentage
import = 0.2

# For each of the import generation select animals within population
# For pop1, import a percentage of males from pop2 and pop3, in pop2 and pop3 mate within population
for (Generation in (1 + nGenerationBurn):(nGenerationEval + nGenerationBurn)) {
  for (pop in c("Pop1", "Pop2", "Pop3")) {

    # If this if the first generation of evaluation, select from burn-in animals
    if (Generation == (1 + nGenerationBurn)) {
      # If this is population 1, create a mixture of pop1, pop2 and pop3 fathers
      if (pop == "Pop1") {
        Females = BaseFemales_pop1
        Males = c(BaseMales_pop1, BaseMales_pop2, BaseMales_pop3)
        
        # If pop == Pop1, mate with a mating plan - the males contains half the import % of males from Pop2 and half the import % of males from Pop3
        # If pop == Pop2, randomly mate
        matingPlan = cbind(rep(Females@id, 2), 
                           c(sample(BaseMales_pop1@id, size = Females@nInd*2*(1-import), replace=T), 
                             sample(BaseMales_pop2@id, size = Females@nInd*2*(import/2), replace=T),
                             sample(BaseMales_pop3@id, size = Females@nInd*2*(import/2), replace=T)
                           ))
        SelCand = makeCross2(females = Females, males = Males, crossPlan = matingPlan, 
                             nProgeny = 2)
      } else if (pop == "Pop2") {
        #If this is Pop2, just select the parents
        Females = BaseFemales_pop2
        Males = BaseMales_pop2
        SelCand = randCross2(females = Females, males = Males, nCrosses = Females@nInd,
                             nProgeny = 2)
      } else if (pop == "Pop3") {
        #If this is Pop2, just select the parents
        Females = BaseFemales_pop3
        Males = BaseMales_pop3
        SelCand = randCross2(females = Females, males = Males, nCrosses = Females@nInd,
                             nProgeny = 2)
      }
    }
    
    #I this is not the first generation of evaluation, select from previous round selected parens
    if (Generation > (1 + nGenerationBurn)) {
      #If this is pop 1, create a mixture of pop1 and pop2 fathers
      if (pop == "Pop1") {
        Females = Females_pop1
        Males = c(Males_pop1, Males_pop2, Males_pop3)
        
        # If pop == Pop1, mate with a mating plan - the males contain the import % of males from Pop2
        # If pop == Pop2, randomly mate
        matingPlan = cbind(rep(Females@id, 2), 
                           c(sample(Males_pop1@id, size = Females@nInd*2*(1-import), replace=T), 
                             sample(c(Males_pop2@id, Males_pop3@id), size = Females@nInd*2*import, replace=T)))
        SelCand = makeCross2(females = Females, males = Males, crossPlan = matingPlan, 
                             nProgeny = 2)
        
      } else if (pop == "Pop2") {
        Females = Females_pop2
        Males = Males_pop2
        SelCand = randCross2(females = Females, males = Males, nCrosses = Females@nInd,
                             nProgeny = 2)
      } else if (pop == "Pop3") {
        Females = Females_pop3
        Males = Males_pop3
        SelCand = randCross2(females = Females, males = Males, nCrosses = Females@nInd,
                             nProgeny = 2)
      }
    }

    # meanG(SelCand)
    
    # Save metrics
    for (gender in c("F", "M")) {
      DataImport[(DataImport$Generation == Generation) & (DataImport$Gender == gender) & (DataImport$Pop == pop), GenMeanCols] =
        colMeans(SelCand@gv[SelCand@gender == gender,])
      DataImport[(DataImport$Generation == Generation) & (DataImport$Gender == gender) & (DataImport$Pop == pop), GenVarCols] =
        diag(var(SelCand@gv[SelCand@gender == gender,]))
    }
    
    # Phenotype
    SelCand = setPheno(pop = SelCand, varE = VarE)

    # Track pedigree
    PedImport = rbind(PedImport,
                    tibble(Generation = Generation,
                           IId        = SelCand@id,
                           FId        = SelCand@father,
                           MId        = SelCand@mother,
                           Gender     = SelCand@gender,
                           Program    = "Eval",
                           Population = pop,
                           PhenoT1    = SelCand@pheno[,1],
                           PhenoT2    = SelCand@pheno[,2],
                           PhenoT3    = SelCand@pheno[,3],
                           TbvT1      = SelCand@gv[, 1],
                           TbvT2      = SelCand@gv[, 2],
                           TbvT3      = SelCand@gv[, 3]))
   
    
    # Compute Accuracies (pheno-TGV) for three populations
    accuracies <- rbind(accuracies, c(pop, Generation, 1, cor(SelCand@pheno[,1], SelCand@gv[,1]) ))
    accuracies <- rbind(accuracies, c(pop, Generation, 2, cor(SelCand@pheno[,2], SelCand@gv[,2]) ))
    accuracies <- rbind(accuracies, c(pop, Generation, 3, cor(SelCand@pheno[,3], SelCand@gv[,3]) ))
    
    # Select within populations
    if (pop == "Pop1") {
      Males_pop1   = selectInd(pop = SelCand, nInd = nMales_pop1,   gender = "M",
                                   use = "pheno", trait = 1)
      Females_pop1 = selectInd(pop = SelCand, nInd = nFemales_pop1, gender = "F")
    } else if (pop == "Pop2") {
      Males_pop2   = selectInd(pop = SelCand, nInd = nMales_pop2,   gender = "M",
                                   use = "pheno", trait = 2)
      Females_pop2 = selectInd(pop = SelCand, nInd = nFemales_pop2, gender = "F")
    } else if (pop == "Pop3") {
      Males_pop3   = selectInd(pop = SelCand, nInd = nMales_pop3,   gender = "M",
                                   use = "pheno", trait = 3)
      Females_pop3 = selectInd(pop = SelCand, nInd = nFemales_pop3, gender = "F")
    }
    # Clean
    rm(SelCand)
  }
}

# Save simulation pedigree
write.table(as.data.frame(PedImport), "SimulatedPedigree.csv", quote=F, row.names=F, sep=",")

# Plot genetic means
rbind(DataBurn, DataImport) %>%
  gather(key = "Metric", value = "Value", GenMeanCols) %>%
  ggplot(., aes(Generation, Value, color = Population)) +
  geom_line() + facet_grid(rows = vars(Metric), cols = vars(Gender)) + 
  ylab(label = "Genetic mean")


# Plot only the one correlated (observed) trait, i.e. trait 1 for pop1, trait2 for pop2, trait 3 for pop3
OneTrait <- rbind(DataBurn, DataImport) %>% gather(key = "Metric", value = "Value", GenMeanCols)
rbind(OneTrait[OneTrait$Population == "Pop1" & OneTrait$Metric == "GenMeanT1",], 
      OneTrait[OneTrait$Population == "Pop2" & OneTrait$Metric == "GenMeanT2",],
      OneTrait[OneTrait$Population == "Pop3" & OneTrait$Metric == "GenMeanT3",]) %>% 
  ggplot(., aes(Generation, Value, color = Metric)) +
  geom_line() + facet_grid(cols = vars(Gender)) + 
  ylab(label = "Genetic mean")

# Plot genetic variances
DataImport %>% gather(key = "Metric", value= "Value", GenVarCols) %>%
  ggplot(., aes(Generation, Value, color = Metric)) +
  geom_line() + facet_grid(rows = vars(Population), cols = vars(Gender)) + 
ylab(label = "Genetic variance")  



##############################################################################
##############################################################################
# ---- Partition analysis ----
###############################################################################
###############################################################################
# Use AlphaPart to partition genetic trends and quantify the sources of genetic gain in pop1

# Create a variable Bv that hold the one correlated (observed) trait, trait 1 for pop1, trait 2 for pop2, trait 3 for pop3
PedImport$Bv <- ifelse(PedImport$Population == "Pop1", PedImport$TbvT1, ifelse(PedImport$Population == "Pop2", PedImport$TbvT2, PedImport$TbvT3))

# Keep only the columns you need for partition analysis
PedImport <- as.data.frame(PedImport[, c("Generation", "IId", "FId", "MId", "Population", "Bv")])

# Rebase the pedigree to hold only the information of the import generations
# This sets generation 11 as the founder population
PedImport <- pedSetBase(PedImport, 
                        PedImport$Generation>10,
                        colId = "IId", colFid = "FId", colMid = "MId")

# We fix or impute missing or erroneous years of birth 
# However, since data is simulated, all the years of birth are known and correct and this is only to demonstrate the function
PedImport <- pedFixBirthYear(PedImport,
                            interval=1,
                            colId = "IId", colFid = "FId", colMid = "MId")

# Partition the breeding values for the observed trait (variable Bv) by population
Part = AlphaPart(x = PedImport,
                 colPath = "Population", 
                 colBV = "Bv",
                 colId = "IId", colFid = "FId", colMid = "MId")

# Summarize the partitioned breeding values by generation
# Keep only the records for population 1
sumPart <- summary(Part, by="Generation", 
                   subset=Part$Bv$Population == "Pop1")

# Plot the summarized results
plot(sumPart, ylab="Mean breeding value")

# We subset the summarized dataset to keep only the contributions of pop1 and pop2
sumPart_subset <- AlphaPartSubset(sumPart, 
                                   paths = c("Pop1", "Pop2"))

# We plot the subset partitions
plot(sumPart_subset, ylab="Mean breeding value")

# We summarise the contributions into the contribution of domestic selection (pop1) and imported gain (pop2 and pop3)
sumPart_sum <- AlphaPartSum(sumPart, 
                            map=list(c("Domestic", "Pop1"), c("Import", "Pop2", "Pop3")))

# Plot the summarised results
plot(sumPart_sum, ylab="Mean breeding value")
