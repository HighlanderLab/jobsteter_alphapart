# ---- Environment ----
library(package = "AlphaSimR")
library(package = "tidyverse")
library(package = "AlphaPart")


# ---- Notes ---- 
### SET THE homeDir to WHERE YOU WANT TO RUN THE SIMULATION
### THE DOWNLOADED Essentials FOLDER SHOULD BE IN THIS DIRECTORY!!!
### For this script to run, you should download renumf90 and blupf90 (Mizstal et al., Univeristy of Georgia) and55+41+ move it into the Essentials folder.
homeDir = getwd()

# ---- Glossary ----
# GN = nucleus
# PN = multiplier
# PN1 = Programme 1 = MaleFlow100
# PN2 = Programme 2 = MaleFlow50

# Base Burn-In (Nucleus)
# Then for each generation: selection in the nucleus and selection in the multiplier
#load the packages

# ---- General simulation parameters ----
nGNMales   =  25
nGNFemales = 500

nPNMales       =  100
nPNFemales     = 4500
pPNFemalesPure = 1/6

nGenerationBurn = 20
nGenerationEval = 20

GenMeanCols = c("GenMeanT1", "GenMeanT2", "GenMeanI")
GenVarCols  = c("GenVarT1",  "GenVarT2",  "GenVarI")

# ---- Base population genomes ----

#simulate a founder population
founderPop = runMacs(nInd = nGNMales + nGNFemales,
                     nChr = 10,
                     segSites = 1000,
                     nThreads = 4,
                     # species = "GENERIC")
                     species = "CATTLE")


# ---- Simulation/Base population parameters ----
# Make folders to hold the information for estimation of breeding values
system('mkdir GN PN1 PN2')
system(paste0('cp ', homeDir, '/Essentials/* ', homeDir, '/PN1/'))
system(paste0('cp ', homeDir, '/Essentials/* ', homeDir, '/PN2/'))
system(paste0('cp ', homeDir, '/Essentials/* ', homeDir, '/GN/'))

# Set a new population
SP = SimParam$new(founderPop)
# Set genetic variance 1 for trait 1 and trait 2
VarA = matrix(data = c(1.0, 0.0, 0.0, 1.0), nrow = 2); cov2cor(VarA)
# Set residual variance 3 for trait 1 and 9 for trait 2, resulting in 0.25 and 0.1 heritability
VarE = matrix(data = c(3.0, 0.0, 0.0, 9.0), nrow = 2); cov2cor(VarE)
# Compute the phenotypic variance
VarP = VarA + VarE; diag(VarA) / diag(VarP)
# Add an aditive trait controled with 1000QTNs per chromosome
SP$addTraitA(nQtlPerChr = 1000, mean = c(0, 0), var = diag(VarA), cor = cov2cor(VarA))
# Set gender randomly
SP$setGender(gender = "yes_rand")
  
# ---- Base Nucleus (GN) population ----
# Set nucleus population
GN = newPop(founderPop)
# Select base nucleus males
BaseGNMales   = GN[GN@gender == "M"]
# Select base nucleus females
BaseGNFemales = GN[GN@gender == "F"]
rm(GN)
  

############ -- BURN-in -- #######################################################################

# Perform nucleus burn-in, select based on index phenotypic value

# Prepare a dataframe to fold the information on generation genetic means
DataBurn = tibble(Generation = rep(1:nGenerationBurn, each=2), Gender = rep(c("F", "M"), nGenerationBurn),
                  GenMeanT1 = NA, GenMeanT2 = NA, GenMeanI = NA,
                  GenVarT1  = NA, GenVarT2  = NA, GenVarI  = NA)
# Prepare a dataframe to hold the pedigree and individuals' information
PedEval <- tibble()

# For each of the 20 generations of burn-in
for (Generation in 1:nGenerationBurn) {

  # Mate --> create selection candidates
  SelCand = randCross2(females = BaseGNFemales, males = BaseGNMales,
                       nCrosses = BaseGNFemales@nInd, nProgeny = 12)
  
  # Save metrics
  for (gender in c("F", "M")) {
    DataBurn[(DataBurn$Generation == Generation) & (DataBurn$Gender == gender), GenMeanCols] =
      c(colMeans(SelCand@gv[SelCand@gender == gender,]),  mean(rowMeans(SelCand@gv[SelCand@gender == gender,])))
    DataBurn[(DataBurn$Generation == Generation) & (DataBurn$Gender == gender), GenVarCols] =
      c(diag(var(SelCand@gv[SelCand@gender == gender,])), var(rowMeans(SelCand@gv[SelCand@gender == gender,])))
  }
  
  
  # Phenotype
  SelCand = setPheno(pop = SelCand, varE = VarE)
  # 
  # if (Generation == 1) {
  #   VarA <- varG(SelCand)
  #   VarE <- varP(SelCand) - varG(SelCand)
  # }
  # 
  # Add the pedigree information on selection candidates
  PedEval = rbind(PedEval,
                  tibble(Generation = Generation,
                         IId        = SelCand@id,
                         FId        = SelCand@father,
                         MId        = SelCand@mother,
                         Gender     = SelCand@gender,
                         Program    = "BurnIn",
                         PhenoT1    = SelCand@pheno[,1],
                         PhenoT2    = SelCand@pheno[,2],
                         EbvT1      = NA,
                         EbvT2      = NA,
                         TbvT1      = SelCand@gv[, 1],
                         TbvT2      = SelCand@gv[, 2]))
  
  # Select nucleus males and females for the new generation
  BaseGNMales   = selectInd(pop = SelCand, nInd = nGNMales,   gender = "M",
                            use = "pheno", trait = function(x) rowMeans(scale(x)))
  BaseGNFemales = selectInd(pop = SelCand, nInd = nGNFemales, gender = "F",
                            use = "pheno", trait = function(x) rowMeans(scale(x)))
}

# Plot genetic means
DataBurn %>%
  gather(key = "Metric", value = "Value", GenMeanCols) %>%
  ggplot(., aes(Generation, Value, color = Metric)) +
  geom_line() +
  ylab(label = "Genetic mean")

# Plot genetic variances
DataBurn %>% gather(key = "Metric", value= "Value", GenVarCols) %>%
  ggplot(., aes(Generation, Value, color = Metric)) +
  geom_line() +
  ylab(label = "Genetic variance")


# ---- Base Multiplier (PN) population ----
# Obtain nucleus individuals
GNInd = c(BaseGNFemales@id, BaseGNMales@id)
# Selection candidates for multiplier females
SelCand = SelCand[!(SelCand@id %in% GNInd)]; SelCand@nInd
# Cheat here - consider all animals are females
SelCand@gender[] = "F"
# Select base multiplier females
BasePNFemales = selectInd(pop = SelCand, nInd = nPNFemales, gender = "F",
                          use = "pheno", trait = function(x) rowMeans(scale(x)))


# Save the burn-in pedigree to use it for both programs
PedEvalBurnIn <- PedEval

# Set a dataframe to hold the nucleus genetic means
DataEvalGN = tibble(Generation = rep(1:nGenerationEval, each=2),
                    Program = NA, Gender = rep(c("M", "F"), nGenerationEval),
                    GenMeanT1 = NA, GenMeanT2 = NA, GenMeanI = NA,
                    GenVarT1  = NA, GenVarT2  = NA, GenVarI  = NA)

# Create data frame for nucleus and both programmes - PN1 (MaleFlow100) and PN2 (MaleFlow50)
DataEvalGN$Program = "GN"
DataEvalPN1 = DataEvalGN
DataEvalPN1$Program = "PN1"
DataEvalPN2 = DataEvalGN
DataEvalPN2$Program = "PN2"

# Create dataframes to hold the information on accuracy of the EBVs in both programmes
accuraciesPN1 <- data.frame(Program = NA, Generation = NA, Trait = NA, Cor = NA)
accuraciesPN2 <- data.frame(Program = NA, Generation = NA, Trait = NA, Cor = NA)

##############################################################################3
# ---- Program 1 (MaleFlow100): PN with 100% GNmales, PN does T1 ----

PedEval <- PedEvalBurnIn  

for (Generation in (1 + nGenerationBurn):(nGenerationEval + nGenerationBurn)) {

  # First select the nucleus
  if (Generation == (1 + nGenerationBurn)) {
    GNFemales = BaseGNFemales
    GNMales = BaseGNMales
    PedEval$Program[PedEval$IId %in% BaseGNFemales@id] <- "GN"
    PedEval$Program[PedEval$IId %in% BaseGNMales@id] <- "GN"
  }
  
  # Mate nucleus animals
  SelCand = randCross2(females = GNFemales, males = GNMales,
                       nCrosses = GNFemales@nInd, nProgeny = 12)
  
  # Save metrics
  for (gender in c("F", "M")) {
    DataEvalGN[(DataEvalGN$Generation == Generation) & (DataEvalGN$Gender == gender), GenMeanCols] =
      c(colMeans(SelCand@gv[SelCand@gender == gender,]),  mean(rowMeans(SelCand@gv[SelCand@gender == gender,])))
    DataEvalGN[(DataEvalGN$Generation == Generation) & (DataEvalGN$Gender == gender), GenVarCols] =
      c(diag(var(SelCand@gv[SelCand@gender == gender,])), var(rowMeans(SelCand@gv[SelCand@gender == gender,])))
  }
  
  # Phenotype
  SelCand = setPheno(pop = SelCand, varE = VarE)

  # Track pedigree for nucleus animals
  PedEval = rbind(PedEval,
                  tibble(Generation = Generation,
                         IId        = SelCand@id,
                         FId        = SelCand@father,
                         MId        = SelCand@mother,
                         Gender     = SelCand@gender,
                         Program    = "GN",
                         PhenoT1    = SelCand@pheno[,1],
                         PhenoT2    = SelCand@pheno[,2],
                         EbvT1      = NA,
                         EbvT2      = NA,
                         TbvT1      = SelCand@gv[, 1],
                         TbvT2      = SelCand@gv[, 2]))
  
  # Estimate EBVs with blupf90
  setwd(paste0(homeDir, "/GN/"))
  print(getwd())
  if (Generation == (1 + nGenerationBurn)) {
    system("rm *ped *dat")
  }
  
  # Create pedigree file
  blupPed <- PedEval[,c("IId", "FId", "MId")]
  blupPed <- blupPed[order(blupPed$IId),]
  write.table(blupPed, "Blupf90.ped", quote=FALSE, row.names=FALSE, col.names=FALSE, sep=" ", na = "0")
  
  # Create phenotype file for trait 1
  blup1dat <- PedEval[PedEval$Program != "BurnIn",c("IId", "PhenoT1")]
  #add a mean column for blupf90
  blup1dat$Mean <- 1
  write.table(blup1dat, "Blupf901.dat",
              quote=FALSE, row.names=FALSE, col.names=FALSE, sep=" ",
              na = "0", append = file.exists("Blupf901.dat"))
  
  
  # Create phenotype file for trait 2
  blup2dat <- PedEval[PedEval$Program == "GN",c("IId", "PhenoT2")]
  #add a mean column for blupf90
  blup2dat$Mean <- 1
  write.table(blup2dat, "Blupf902.dat",
              quote=FALSE, row.names=FALSE, col.names=FALSE, sep=" ",
              na = "0", append = file.exists("Blupf902.dat"))
  
  # Set residual values in renumf90
  system(paste0('sed "s/ResVariance/', diag(VarE)[1], '/g" renumf901_generic.par > renumf901.par '))
  system(paste0('sed "s/ResVariance/', diag(VarE)[2], '/g" renumf902_generic.par > renumf902.par '))
  system(paste0('sed -i "s/GenVariance/', diag(VarA)[1], '/g" renumf901.par '))
  system(paste0('sed -i "s/GenVariance/', diag(VarA)[2], '/g" renumf902.par '))
  
  # Evaluate Trait 1
  # Run the program to renum the pedigree and prepare files (download renumf90 from Ignacy Misztal - University of Georgia site)
  system(command = "./renumf90 < renumParam1")
  # Run evalution (download blupf90 from Ignacy Misztal - University of Georgia site)
  system(command = "./blupf90 renf90.par")
  # Rename the IDs back to the original
  system(command = "bash Match_AFTERRenum.sh")
  # Rename the solutions file to include the trait (1) and generation number
  system(command = paste0("mv renumbered_Solutions renumbered_Solutions_1_", Generation))
  
  sol1 <- read.table(paste0("renumbered_Solutions_1_", Generation))[,2:3]
  sol1 <- sol1[sol1$V2 %in% PedEval$IId,]
  sol1 <- sol1[order(match(sol1$V2, PedEval$IId)),]
  PedEval$EbvT1 <- sol1$V3
  
  # Evaluate Trait 2 - same as for trait 1
  system(command = "./renumf90 < renumParam2")
  system(command = "./blupf90 renf90.par")
  system(command = "bash Match_AFTERRenum.sh")
  system(command = paste0("mv renumbered_Solutions renumbered_Solutions_2_", Generation))
  
  sol2 <- read.table(paste0("renumbered_Solutions_2_", Generation))[,2:3]
  sol2 <- sol2[sol2$V2 %in% PedEval$IId,]
  sol2 <- sol2[order(match(sol2$V2, PedEval$IId)),]
  PedEval$EbvT2 <- sol2$V3
  
  
  # Set EBVs for SelCand
  SelCand@ebv = as.matrix(PedEval[PedEval$IId %in% SelCand@id, c("EbvT1", "EbvT2")])
  accuraciesPN1 <- rbind(accuraciesPN1, c("GN", Generation, 1, cor(SelCand@ebv[,1], SelCand@gv[,1]) ))
  accuraciesPN1 <- rbind(accuraciesPN1, c("GN", Generation, 2, cor(SelCand@ebv[,2], SelCand@gv[,2]) ))
  
  # Select
  print(summary(PedEval$EbvT1))
  print(summary(PedEval$EbvT2))
  print("Selecting GN1 animals")
  GNMales   = selectInd(pop = SelCand, nInd = nGNMales,   gender = "M",
                        use = "ebv", trait = function(x) rowMeans(scale(x)))
  GNFemales = selectInd(pop = SelCand, nInd = nGNFemales, gender = "F",
                        use = "ebv", trait = function(x) rowMeans(scale(x)))
  # Clean
  rm(SelCand)
  

  # Next select in the multiplier
  if (Generation == (nGenerationBurn + 1)) {
    # If this is the first generation, the mother are the base multiplier females
    # Otherwise they are the selected multiplier females from previous generation
    PNFemales1 = BasePNFemales
    PedEval$Program[PedEval$IId %in% BasePNFemales@id] <- "PN1"
    PNFemales1ForPure = selectInd(pop = PNFemales1, nInd = PNFemales1@nInd * pPNFemalesPure,
                                  use = "pheno", trait = function(x) rowMeans(scale(x), na.rm = TRUE))
  } else {
    PNFemales1ForPure = selectInd(pop = PNFemales1, nInd = PNFemales1@nInd * pPNFemalesPure,
                                  use = "ebv", trait = function(x) rowMeans(scale(x), na.rm = TRUE))
  }
  
  # Mate Multiplier females and nucleus males
  SelCand = randCross2(females = PNFemales1ForPure, males = GNMales,
                       nCrosses = PNFemales1ForPure@nInd, nProgeny = 12)

  
  # Save metrics
  for (gender in c("F", "M")) {
    DataEvalPN1[(DataEvalPN1$Generation == Generation) & (DataEvalPN1$Gender == gender), GenMeanCols] =
      c(colMeans(SelCand@gv[SelCand@gender == gender,]),  mean(rowMeans(SelCand@gv[SelCand@gender == gender,])))
    DataEvalPN1[(DataEvalPN1$Generation == Generation) & (DataEvalPN1$Gender == gender), GenVarCols] =
      c(diag(var(SelCand@gv[SelCand@gender == gender,])), var(rowMeans(SelCand@gv[SelCand@gender == gender,])))
  }
  
  
  # Phenotype, remove the phenotype for trait 2 (not measured in the multiplier)
  SelCand = setPheno(pop = SelCand, varE = VarE)
  SelCand@pheno[, 2] <- NA
  
  
  # Track pedigree, multiplier is +1 generation
  PedEval = rbind(PedEval, 
                  tibble(Generation = Generation +1,
                         IId        = SelCand@id,
                         FId        = SelCand@father,
                         MId        = SelCand@mother,
                         Gender     = SelCand@gender,
                         Program    = "PN1",
                         PhenoT1    = SelCand@pheno[,1],
                         PhenoT2    = SelCand@pheno[,2],
                         EbvT1      = NA,
                         EbvT2      = NA,
                         TbvT1      = SelCand@gv[, 1],
                         TbvT2      = SelCand@gv[, 2]))
  
  
  # Estimate EBVs with blupf90
  setwd(paste0(homeDir, "/PN1/"))

  if (Generation == (1 + nGenerationBurn)) {
    system("rm *ped *dat")
  }
  
  # Create pedigree file
  blupPed <- PedEval[,c("IId", "FId", "MId")]
  blupPed <- blupPed[order(blupPed$IId),]
  write.table(blupPed, "Blupf90.ped", quote=FALSE, row.names=FALSE, col.names=FALSE, sep=" ", na = "0")
  
  # Create phenotype file for trait 1
  blup1dat <- PedEval[PedEval$Program != "BurnIn",c("IId", "PhenoT1")]
  # Add a mean column for blupf90
  blup1dat$Mean <- 1
  write.table(blup1dat, "Blupf901.dat",
              quote=FALSE, row.names=FALSE, col.names=FALSE, sep=" ",
              na = "0", append = file.exists("Blupf901.dat"))
  
  # Create phenotype file for trait 2
  blup2dat <- PedEval[PedEval$Program == "GN",c("IId", "PhenoT2")]
  blup2dat$Mean <- 1
  write.table(blup2dat, "Blupf902.dat",
              quote=FALSE, row.names=FALSE, col.names=FALSE, sep=" ",
              na = "0", append = file.exists("Blupf902.dat"))   

  # Set residual values in renumf90
  system(paste0('sed "s/ResVariance/', diag(VarE)[1], '/g" renumf901_generic.par > renumf901.par '))
  system(paste0('sed "s/ResVariance/', diag(VarE)[2], '/g" renumf902_generic.par > renumf902.par '))
  system(paste0('sed -i "s/GenVariance/', diag(VarA)[1], '/g" renumf901.par '))
  system(paste0('sed -i "s/GenVariance/', diag(VarA)[2], '/g" renumf902.par '))
  
  # Evaluate Trait 1
  # Run the program to renum the pedigree and prepare files
  system(command = "./renumf90 < renumParam1")
  # Run evalution
  system(command = "./blupf90 renf90.par")
  # Rename the IDs back to the original
  system(command = "bash Match_AFTERRenum.sh")
  # Rename the solutions file to include the trait (1) and generation number
  system(command = paste0("mv renumbered_Solutions renumbered_Solutions_1_", Generation))
  
  # Read in the solutions and bind them to the pedigree
  sol1 <- read.table(paste0("renumbered_Solutions_1_", Generation))[,2:3]
  sol1 <- sol1[sol1$V2 %in% PedEval$IId, ]
  sol1 <- sol1[order(match(sol1$V2, PedEval$IId)),]
  PedEval$EbvT1 <- sol1$V3  
  
  # Evaluate Trait 2 - same as for trait 1
  system(command = "./renumf90 < renumParam2")
  system(command = "./blupf90 renf90.par")
  system(command = "bash Match_AFTERRenum.sh")
  system(command = paste0("mv renumbered_Solutions renumbered_Solutions_2_", Generation))
  
  sol2 <- read.table(paste0("renumbered_Solutions_2_", Generation))[,2:3]
  sol2 <- sol2[sol2$V2 %in% PedEval$IId, ]
  sol2 <- sol2[order(match(sol2$V2, PedEval$IId)),]
  PedEval$EbvT2 <- sol2$V3

  # Set EBVs for SelCand
  SelCand@ebv = as.matrix(PedEval[PedEval$IId %in% SelCand@id, c("EbvT1", "EbvT2")])
  # Track the accuracies
  accuraciesPN1 <- rbind(accuraciesPN1, c("PN1", Generation, 1, cor(SelCand@ebv[,1], SelCand@gv[,1]) ))
  accuraciesPN1 <- rbind(accuraciesPN1, c("PN1", Generation, 2, cor(SelCand@ebv[,2], SelCand@gv[,2]) ))
  
  # Select multiplier females and males (although you do not use males in this programme)
  # PNMales1   = selectInd(pop = SelCand, nInd = nPNMales,   gender = "M",
  #                      use = "ebv", trait = function(x) rowMeans(x, na.rm = TRUE))
  PNFemales1 = selectInd(pop = SelCand, nInd = nPNFemales, gender = "F",
                        use = "ebv", trait = function(x) rowMeans(scale(x), na.rm = TRUE))
  
  # Clean
  rm(SelCand)
}

# Assign the tracked pedigree as PedEval1, same for DataEval
PedEval1 = PedEval
rm(PedEval)
DataEvalGN1 = DataEvalGN
rm(DataEvalGN)

# Write the files
write.table(PedEval1, paste0("PedEval1.csv"), quote=FALSE, row.names=FALSE)
write.table(DataEvalGN1, paste0("DataEval1.csv"), quote=FALSE, row.names=FALSE)
write.table(accuraciesPN1, paste0("Accuracies1.csv"), quote=FALSE, row.names=FALSE)



##############################################################################
##############################################################################
# ---- Program 2 = MaleFlow50: PN with 50% GNmales, PN does T1 ----

# Set the DataFrame to track the nucleus information
DataEvalGN = tibble(Generation = rep(1:nGenerationEval, each=2),
                    Program = NA, Gender = rep(c("M", "F"), nGenerationEval),
                    GenMeanT1 = NA, GenMeanT2 = NA, GenMeanI = NA,
                    GenVarT1  = NA, GenVarT2  = NA, GenVarI  = NA)
DataEvalGN$Program = "GN"

# Obtain the burn-in pedigree
PedEval <- PedEvalBurnIn

# Select from generation 20 to 40
for (Generation in (1 + nGenerationBurn):(nGenerationEval + nGenerationBurn)) {
  # If this is the first generation, the parents are the base animals
  if (Generation == (1 + nGenerationBurn)) {
    GNFemales = BaseGNFemales
    GNMales = BaseGNMales
    PedEval$Program[PedEval$IId %in% BaseGNFemales@id] <- "GN"
    PedEval$Program[PedEval$IId %in% BaseGNMales@id] <- "GN"
  }
  
  # Mate nucleus animals
  SelCand = randCross2(females = GNFemales, males = GNMales,
                       nCrosses = GNFemales@nInd, nProgeny = 12)
  
  # Save metrics
  for (gender in c("F", "M")) {
    DataEvalGN[(DataEvalGN$Generation == Generation) & (DataEvalGN$Gender == gender), GenMeanCols] =
      c(colMeans(SelCand@gv[SelCand@gender == gender,]),  mean(rowMeans(SelCand@gv[SelCand@gender == gender,])))
    DataEvalGN[(DataEvalGN$Generation == Generation) & (DataEvalGN$Gender == gender), GenVarCols] =
      c(diag(var(SelCand@gv[SelCand@gender == gender,])), var(rowMeans(SelCand@gv[SelCand@gender == gender,])))
  }
  
  # Phenotype
  SelCand = setPheno(pop = SelCand, varE = VarE)
  
  # Track pedigree
  PedEval = rbind(PedEval,
                  tibble(Generation = Generation,
                         IId        = SelCand@id,
                         FId        = SelCand@father,
                         MId        = SelCand@mother,
                         Gender     = SelCand@gender,
                         Program    = "GN",
                         PhenoT1    = SelCand@pheno[,1],
                         PhenoT2    = SelCand@pheno[,2],
                         EbvT1      = NA,
                         EbvT2      = NA,
                         TbvT1      = SelCand@gv[, 1],
                         TbvT2      = SelCand@gv[, 2]))
  
  
  # Estimate EBVs with blupf90
  setwd(paste0(homeDir, "/GN/"))
  
  if (Generation == (1 + nGenerationBurn)) {
    system("rm *ped *dat")
  }
  
  # Create pedigree file
  blupPed <- PedEval[,c("IId", "FId", "MId")]
  blupPed <- blupPed[order(blupPed$IId),]
  write.table(blupPed, "Blupf90.ped", quote=FALSE, row.names=FALSE, col.names=FALSE, sep=" ", na = "0")
  
  # Create phenotype file for trait 1
  # Trait 1
  blup1dat <- PedEval[PedEval$Program != "BurnIn", c("IId", "PhenoT1")]
  # Add a mean column for blupf90
  blup1dat$Mean <- 1
  write.table(blup1dat, "Blupf901.dat",
              quote=FALSE, row.names=FALSE, col.names=FALSE, sep=" ",
              na = "0", append = file.exists("Blupf901.dat"))


  # Create phenotype file for trait 2
  blup2dat <- PedEval[PedEval$Program == "GN",c("IId", "PhenoT2")]
  blup2dat$Mean <- 1
  write.table(blup2dat, "Blupf902.dat",
              quote=FALSE, row.names=FALSE, col.names=FALSE, sep=" ",
              na = "0", append = file.exists("Blupf902.dat"))    

  ## set residual values in renumf90
  system(paste0('sed "s/ResVariance/', diag(VarE)[1], '/g" renumf901_generic.par > renumf901.par '))
  system(paste0('sed "s/ResVariance/', diag(VarE)[2], '/g" renumf902_generic.par > renumf902.par '))
  system(paste0('sed -i "s/GenVariance/', diag(VarA)[1], '/g" renumf901.par '))
  system(paste0('sed -i "s/GenVariance/', diag(VarA)[2], '/g" renumf902.par '))
  
  # Evaluate Trait 1
  # Run the program to renum the pedigree and prepare files
  system(command = "./renumf90 < renumParam1")
  # Run evalution
  system(command = "./blupf90 renf90.par")
  # Rename the IDs back to the original
  system(command = "bash Match_AFTERRenum.sh")
  # Rename the solutions file to include the trait (1) and generation number
  system(command = paste0("mv renumbered_Solutions renumbered_Solutions_1_", Generation))
  
  # Read in the solutions and bind them to the pedigree
  sol1 <- read.table(paste0("renumbered_Solutions_1_", Generation))[,2:3]
  sol1 <- sol1[sol1$V2 %in% PedEval$IId, ]
  sol1 <- sol1[order(match(sol1$V2, PedEval$IId)),]
  PedEval$EbvT1 <- sol1$V3
  
  # Evaluate Trait 2 - same as for trait 1
  system(command = "./renumf90 < renumParam2")
  system(command = "./blupf90 renf90.par")
  system(command = "bash Match_AFTERRenum.sh")
  system(command = paste0("mv renumbered_Solutions renumbered_Solutions_2_", Generation))
  
  sol2 <- read.table(paste0("renumbered_Solutions_2_", Generation))[,2:3]
  sol2 <- sol2[sol2$V2 %in% PedEval$IId, ]
  sol2 <- sol2[order(match(sol2$V2, PedEval$IId)),]
  PedEval$EbvT2 <- sol2$V3
  
  
  # Set EBVs for SelCand
  SelCand@ebv = as.matrix(PedEval[PedEval$IId %in% SelCand@id, c("EbvT1", "EbvT2")])
  accuraciesPN2 <- rbind(accuraciesPN2, c("GN", Generation, 1, cor(SelCand@ebv[,1], SelCand@gv[,1]) ))
  accuraciesPN2 <- rbind(accuraciesPN2, c("GN", Generation, 2, cor(SelCand@ebv[,2], SelCand@gv[,2]) ))
  
  # Select
  GNMales   = selectInd(pop = SelCand, nInd = nGNMales,   gender = "M",
                        use = "ebv", trait = function(x) rowMeans(scale(x)))
  GNFemales = selectInd(pop = SelCand, nInd = nGNFemales, gender = "F",
                        use = "ebv", trait = function(x) rowMeans(scale(x)))
  # Clean
  rm(SelCand)
 
  
  # Next select in the multiplier
  if (Generation == (1 + nGenerationBurn)) {
    # If this is the first generation of selection, the mother are the base multiplier females
    # Otherwise they are the selected multiplier females from previous generation
    PNFemales2 = BasePNFemales
    PNMales2 = BaseGNMales
    PedEval$Program[PedEval$IId %in% BasePNFemales@id] <- "PN2"
    PNFemales2ForPure = selectInd(pop = PNFemales2, nInd = PNFemales2@nInd * pPNFemalesPure,
                                  use = "pheno", trait = function(x) rowMeans(scale(x), na.rm = TRUE))
  } else {
    PNFemales2ForPure = selectInd(pop = PNFemales2, nInd = PNFemales2@nInd * pPNFemalesPure,
                                use = "ebv", trait = function(x) rowMeans(scale(x), na.rm = TRUE))
  }
  SelCand = randCross2(females = PNFemales2ForPure, males = c(GNMales, PNMales2),
                       nCrosses = PNFemales2ForPure@nInd, nProgeny = 12)

  # Save metrics
  for (gender in c("F", "M")) {
    DataEvalPN2[(DataEvalPN2$Generation == Generation) & (DataEvalPN2$Gender == gender), GenMeanCols] =
      c(colMeans(SelCand@gv[SelCand@gender == gender,]),  mean(rowMeans(SelCand@gv[SelCand@gender == gender,])))
    DataEvalPN2[(DataEvalPN2$Generation == Generation) & (DataEvalPN2$Gender == gender), GenVarCols] =
      c(diag(var(SelCand@gv[SelCand@gender == gender,])), var(rowMeans(SelCand@gv[SelCand@gender == gender,])))
  }
  
  # Phenotype, remove the phenotype for trait 2 (not measured in the multiplier)
  SelCand = setPheno(pop = SelCand, varE = VarE)
  SelCand@pheno[, 2] = NA
  
  
  # Track pedigree, multiplier is +1 generation
  PedEval = rbind(PedEval, 
                  tibble(Generation = Generation + 1,
                         IId        = SelCand@id,
                         FId        = SelCand@father,
                         MId        = SelCand@mother,
                         Gender     = SelCand@gender,
                         Program    = "PN2",
                         PhenoT1    = SelCand@pheno[,1],
                         PhenoT2    = SelCand@pheno[,2],
                         EbvT1      = NA,
                         EbvT2      = NA,
                         TbvT1      = SelCand@gv[, 1],
                         TbvT2      = SelCand@gv[, 2]))
  
  # Estimate EBVs with blupf90
  setwd(paste0(homeDir, "/PN2/"))
  if (Generation == (1 + nGenerationBurn)) {
    system("rm *ped *dat")
  }
  
  # Create pedigree file
  blupPed <- PedEval[,c("IId", "FId", "MId")]
  blupPed$FId[is.na(blupPed$FId)] <- 0
  blupPed$MId[is.na(blupPed$MId)] <- 0
  blupPed <- blupPed[order(blupPed$IId),]
  write.table(blupPed, "Blupf90.ped", quote=FALSE, row.names=FALSE, col.names=FALSE, sep=" ", na = "0")
  
  # Create phenotype file for trait 1
  blup1dat <- PedEval[PedEval$Program != "BurnIn",c("IId", "PhenoT1")]
  # Add a mean column for blupf90
  blup1dat$Mean <- 1
  write.table(blup1dat, "Blupf901.dat",
              quote=FALSE, row.names=FALSE, col.names=FALSE, sep=" ",
              na = "0", append = file.exists("Blupf901.dat"))
  
  
  # Create phenotype file for trait 2
  blup2dat <- PedEval[PedEval$Program == "GN",c("IId", "PhenoT2")]
  blup2dat$Mean <- 1
  write.table(blup2dat, "Blupf902.dat",
              quote=FALSE, row.names=FALSE, col.names=FALSE, sep=" ",
              na = "0", append = file.exists("Blupf902.dat"))    


  # Set residual values in renumf90
  system(paste0('sed "s/ResVariance/', diag(VarE)[1], '/g" renumf901_generic.par > renumf901.par '))
  system(paste0('sed "s/ResVariance/', diag(VarE)[2], '/g" renumf902_generic.par > renumf902.par '))
  system(paste0('sed -i "s/GenVariance/', diag(VarA)[1], '/g" renumf901.par '))
  system(paste0('sed -i "s/GenVariance/', diag(VarA)[2], '/g" renumf902.par '))
  
  # Evaluate Trait 1
  # Run the program to renum the pedigree and prepare files
  system(command = "./renumf90 < renumParam1")
  # Run evalution
  system(command = "./blupf90 renf90.par")
  system(command = "./blupf90 renf90.par")
  # Rename the IDs back to the original
  system(command = "bash Match_AFTERRenum.sh")
  # Rename the solutions file to include the trait (1) and generation number
  system(command = paste0("mv renumbered_Solutions renumbered_Solutions_1_", Generation))
  
  # Read in the solutions and bind them to the pedigree
  sol1 <- read.table(paste0("renumbered_Solutions_1_", Generation))[,2:3]
  sol1 <- sol1[sol1$V2 %in% PedEval$IId, ]
  sol1 <- sol1[order(match(sol1$V2, PedEval$IId)),]
  PedEval$EbvT1 <- sol1$V3
  
  # Evaluate Trait 2 - same as for trait 1
  system(command = "./renumf90 < renumParam2")
  system(command = "./blupf90 renf90.par")
  system(command = "bash Match_AFTERRenum.sh")
  system(command = paste0("mv renumbered_Solutions renumbered_Solutions_2_", Generation))
  
  sol2 <- read.table(paste0("renumbered_Solutions_2_", Generation))[,2:3]
  sol2 <- sol2[sol2$V2 %in% PedEval$IId, ]
  sol2 <- sol2[order(match(sol2$V2, PedEval$IId)),]
  PedEval$EbvT2 <- sol2$V3
  
  # Set EBVs for SelCand
  SelCand@ebv = as.matrix(PedEval[PedEval$IId %in% SelCand@id, c("EbvT1", "EbvT2")])
  accuraciesPN2 <- rbind(accuraciesPN2, c("PN2", Generation, 1, cor(SelCand@ebv[,1], SelCand@gv[,1]) ))
  accuraciesPN2 <- rbind(accuraciesPN2, c("PN2", Generation, 2, cor(SelCand@ebv[,2], SelCand@gv[,2]) ))
  
  # Select
  PNMales2   = selectInd(pop = SelCand, nInd = nPNMales,   gender = "M",
                         use = "ebv", trait = function(x) rowMeans(scale(x), na.rm = TRUE))
  PNFemales2 = selectInd(pop = SelCand, nInd = nPNFemales, gender = "F",
                         use = "ebv", trait = function(x) rowMeans(scale(x), na.rm = TRUE))
  

  # Clean
  rm(SelCand)
}

# Assign the tracked pedigree as PedEval2, same for DataEval
PedEval2 = PedEval
rm(PedEval)
DataEvalGN2 <- DataEvalGN
rm(DataEvalGN)

# Write the files
accuraciesPN2$Cor <- as.numeric(accuraciesPN2$Cor)
write.table(PedEval2, paste0("PedEval2.csv"), quote=FALSE, row.names=FALSE)
write.table(DataEvalGN2, paste0("DataEval2.csv"), quote=FALSE, row.names=FALSE)
write.table(accuraciesPN2, paste0("Accuracies2.csv"), quote=FALSE, row.names=FALSE)



