library(AlphaPart)

#set the homeDir to wherever your simulation results are
homeDir = ""

# Programme 1 - MaleFlow100
PedEval1 <- read.table(paste0(homeDir, "/PN1/PedEval1.csv"), header=TRUE)
PedEval1 <- PedEval1[PedEval1$Generation %in% 20:41,]
#Burn-in was only nucleus
PedEval1$Program[PedEval1$Program == "BurnIn"] <- "GN"

# Programme 2 - MaleFlow50
PedEval2 <- read.table(paste0(homeDir, "/PedEval2.csv"), header=TRUE)
PedEval2 <- PedEval2[PedEval2$Generation %in% 20:41,]
PedEval2$Program[PedEval2$Program == "BurnIn"] <- "GN"


###################################################################
# 1) PARTITION THE COMPLETE PEDIGREE
###################################################################
# 1. 1) Partition Estimated Breeding Values (EBVs)
  
# ---- PN1 ----
# Compute the index EBV
PedEval1$EbvI = 0.5 * (PedEval1$EbvT1 + PedEval1$EbvT2)
# Standardise onto generation 20
PedEval1$EbvT1_s <- (PedEval1$EbvT1 - mean(PedEval1$EbvT1[PedEval1$Generation == 20])) / sd(PedEval1$TbvT1[PedEval1$Generation == 20])
PedEval1$EbvT2_s <- (PedEval1$EbvT2 - mean(PedEval1$EbvT2[PedEval1$Generation == 20])) / sd(PedEval1$TbvT2[PedEval1$Generation == 20])
# Compute standardised index
PedEval1$EbvI_s <- 0.5 * (PedEval1$EbvT1_s + PedEval1$EbvT2_s)
  
# ---- PN2 ----
# Compute the index EBV
PedEval2$EbvI = 0.5 * (PedEval2$EbvT1 + PedEval2$EbvT2)
# Standardise onto generation 20
PedEval2$EbvT1_s <- (PedEval2$EbvT1 - mean(PedEval2$EbvT1[PedEval2$Generation == 20])) / sd(PedEval2$TbvT1[PedEval2$Generation == 20])
PedEval2$EbvT2_s <- (PedEval2$EbvT2 - mean(PedEval2$EbvT2[PedEval2$Generation == 20])) / sd(PedEval2$TbvT2[PedEval2$Generation == 20])
# Compute standardised index
PedEval2$EbvI_s <- 0.5 * (PedEval2$EbvT1_s + PedEval2$EbvT2_s)
  
#Create Tier-Gender variable
PedEval1$ProgramGender = paste(PedEval1$Program, PedEval1$Gender, sep = "-")
PedEval2$ProgramGender = paste(PedEval2$Program, PedEval2$Gender, sep = "-")

# Partition Programme 1 - MaleFlow100
Part1 = AlphaPart(x = as.data.frame(PedEval1), sort = FALSE,
                  colId = "IId", colFid = "FId", colMid = "MId",
                  colPath = "ProgramGender", colAGV = c("EbvT1_s", "EbvT2_s", "EbvI_s"))

# Partition Programme 2 - MaleFlow50
Part2 = AlphaPart(x = as.data.frame(PedEval2), sort = FALSE,
                  colId = "IId", colFid = "FId", colMid = "MId",
                  colPath = "ProgramGender", colAGV = c("EbvT1_s", "EbvT2_s", "EbvI_s"))
 
# Plot the results
plot(summary(Part1, by="Generation"))
plot(summary(Part2, by="Generation"))

###################################################################
# 1.2) Partition True Breeding Values (TBVs)

# ---- PN1 ----
# Compute the index TBV
PedEval1$TbvI = 0.5 * (PedEval1$TbvT1 + PedEval1$TbvT2)
# Standardise onto generation 20
PedEval1$TbvT1_s <- (PedEval1$TbvT1 -  mean(PedEval1$TbvT1[PedEval1$Generation == 20])) / sd(PedEval1$TbvT1[PedEval1$Generation == 20])
PedEval1$TbvT2_s <- (PedEval1$TbvT2 -  mean(PedEval1$TbvT2[PedEval1$Generation == 20])) / sd(PedEval1$TbvT2[PedEval1$Generation == 20])
# Compute standardised index
PedEval1$TbvI_s <- 0.5 * (PedEval1$TbvT1_s + PedEval1$TbvT2_s)


# ---- PN2 ----
# Compute the index TBV
PedEval2$TbvI = 0.5 * (PedEval2$TbvT1 + PedEval2$TbvT2)
# Standardise onto generation 20
PedEval2$TbvT1_s <- (PedEval2$TbvT1 -  mean(PedEval2$TbvT1[PedEval1$Generation == 20])) / sd(PedEval2$TbvT1[PedEval1$Generation == 20])
PedEval2$TbvT2_s <- (PedEval2$TbvT2 -  mean(PedEval2$TbvT2[PedEval1$Generation == 20])) / sd(PedEval2$TbvT2[PedEval1$Generation == 20])
# Compute standardised index
PedEval2$TbvI_s <- 0.5 * (PedEval2$TbvT1_s + PedEval2$TbvT2_s)

# Create tier-gender variable
PedEval1$ProgramGender = paste(PedEval1$Program, PedEval1$Gender, sep = "-")
PedEval2$ProgramGender = paste(PedEval2$Program, PedEval2$Gender, sep = "-")

# Partition Programme 1 - MaleFlow100
Part1g = AlphaPart(x = as.data.frame(PedEval1), sort = FALSE,
                  colId = "IId", colFid = "FId", colMid = "MId",
                  colPath = "ProgramGender", colAGV = c("TbvT1_s", "TbvT2_s", "TbvI_s"))

# Partition Programme 2 - MaleFlow50
Part2g = AlphaPart(x = as.data.frame(PedEval2), sort = FALSE,
                  colId = "IId", colFid = "FId", colMid = "MId",
                  colPath = "ProgramGender", colAGV = c("TbvT1_s", "TbvT2_s", "TbvI_s"))
  
# Plot the results
plot(summary(Part1g, by="Generation"))
plot(summary(Part2g, by="Generation"))

###################################################################
# 2) PARTITION NUCLEUS / MULTIPLIER GENETIC TREND SEPARATELY
###################################################################
# 2. 1) Partition Estimated Breeding Values (EBVs)

# Programme 1
# Extract nucleus partition
Part1GN <- Part1
Part1GN$EbvT1_s <- Part1GN$EbvT1_s[Part1GN$EbvT1_s$Program == "GN",]
Part1GN$EbvT2_s <- Part1GN$EbvT2_s[Part1GN$EbvT2_s$Program == "GN",]
Part1GN$EbvI_s <- Part1GN$EbvI_s[Part1GN$EbvI_s$Program == "GN",]

# Extract multiplier partition
Part1PN <- Part1
Part1PN$EbvT1_s <- Part1PN$EbvT1_s[Part1PN$EbvT1_s$Program == "PN1",]
Part1PN$EbvT2_s <- Part1PN$EbvT2_s[Part1PN$EbvT2_s$Program == "PN1",]
Part1PN$EbvI_s <- Part1PN$EbvI_s[Part1PN$EbvI_s$Program == "PN1",]

# Programme 2
# Extract nucleus partition
Part2GN <- Part2
Part2GN$EbvT1_s <- Part2GN$EbvT1_s[Part2GN$EbvT1_s$Program == "GN",]
Part2GN$EbvT2_s <- Part2GN$EbvT2_s[Part2GN$EbvT2_s$Program == "GN",]
Part2GN$EbvI_s <- Part2GN$EbvI_s[Part2GN$EbvI_s$Program == "GN",]

# Extract multiplier partition
Part2PN <- Part2
Part2PN$EbvT1_s <- Part2PN$EbvT1_s[Part2PN$EbvT1_s$Program == "PN2",]
Part2PN$EbvT2_s <- Part2PN$EbvT2_s[Part2PN$EbvT2_s$Program == "PN2",]
Part2PN$EbvI_s <- Part2PN$EbvI_s[Part2PN$EbvI_s$Program == "PN2",]


###################################################################
# 2. 2) Partition True Breeding Values (TBVs)

# Programme 1
# Extract nucleus partition
Part1gGN <- Part1g
Part1gGN$TbvT1_s <- Part1gGN$TbvT1_s[Part1gGN$TbvT1_s$Program == "GN",]
Part1gGN$TbvT2_s <- Part1gGN$TbvT2_s[Part1gGN$TbvT2_s$Program == "GN",]
Part1gGN$TbvI_s <- Part1gGN$TbvI_s[Part1gGN$TbvI_s$Program == "GN",]


# Extract multiplier partition
Part1gPN <- Part1g
Part1gPN$TbvT1_s <- Part1gPN$TbvT1_s[Part1gPN$TbvT1_s$Program == "PN1",]
Part1gPN$TbvT2_s <- Part1gPN$TbvT2_s[Part1gPN$TbvT2_s$Program == "PN1",]
Part1gPN$TbvI_s <- Part1gPN$TbvI_s[Part1gPN$TbvI_s$Program == "PN1",]

# Programme 2
# Extract nucleus partition
Part2gGN <- Part2g
Part2gGN$TbvT1_s <- Part2gGN$TbvT1_s[Part2gGN$TbvT1_s$Program == "GN",]
Part2gGN$TbvT2_s <- Part2gGN$TbvT2_s[Part2gGN$TbvT2_s$Program == "GN",]
Part2gGN$TbvI_s <- Part2gGN$TbvI_s[Part2gGN$TbvI_s$Program == "GN",]

# Extract multiplier partition
Part2gPN <- Part2g
Part2gPN$TbvT1_s <- Part2gPN$TbvT1_s[Part2gPN$TbvT1_s$Program == "PN2",]
Part2gPN$TbvT2_s <- Part2gPN$TbvT2_s[Part2gPN$TbvT2_s$Program == "PN2",]
Part2gPN$TbvI_s <- Part2gPN$TbvI_s[Part2gPN$TbvI_s$Program == "PN2",]


###################################################################
# 2. 3) Summarise EBV partitions by Generation to obtain partial genetic trend for the paths (tier-gender)

# Programme 1, EBV
# Nucleus
Part1GNSummary = summary(object = Part1GN, by = "Generation")
# Plot the results
plot(Part1GNSummary)
# Multiplier
Part1PNSummary = summary(object = Part1PN, by = "Generation")
# Plot the results
plot(Part1PNSummary)

# Programme 2, EBV
# Nucleus
Part2GNSummary = summary(object = Part2GN, by = "Generation")
# Plot the results
plot(Part2GNSummary)
# Multiplier
Part2PNSummary = summary(object = Part2PN, by = "Generation")
# Plot the results
plot(Part2PNSummary)

###################################################################
# 2. 4) Summarise TBV partitions by Generation to obtain partial genetic trend for the paths (tier-gender)
# Programme 1, TBV
# Nucleus
Part1gGNSummary = summary(object = Part1gGN, by = "Generation")
# Plot the results
plot(Part1gGNSummary)
# Multiplier
Part1gPNSummary = summary(object = Part1gPN, by = "Generation")
# Plot the results
plot(Part1gPNSummary)

# Programme 2, TBV
# Nucleus
Part2gGNSummary = summary(object = Part2gGN, by = "Generation")
# Plot the results
plot(Part2gGNSummary)
# Multiplier
Part2gPNSummary = summary(object = Part2gPN, by = "Generation")
# Plot the results
plot(Part2gPNSummary)
  
###################################################################
# 2. 5) Bind partitions for the two programmes, Ebv / Tbv, and two trait + index
# EbvT1, EbvT2 and EbvI
# TbvT1, TbvT2 and TbvI for TBV
partitionDF1 = data.frame()
partitionDF2 = data.frame()

for (trait in c("T1", "T2", "I")) {
  # NUCLEUS
  # Programme 1
  # EBV
  t1 <- Part1GNSummary[[paste0("Ebv", trait, "_s")]]$abs
  t1$Program <- "PN1"
  t1$Trait <- trait
  t1$value <- "Ebv"
  t1$Population <- "GN1"
  partitionDF1 <- rbind(partitionDF1, t1)
  
  # TBV
  t1g <- Part1gGNSummary[[paste0("Tbv", trait, "_s")]]$abs
  t1g$Program <- "PN1"
  t1g$Trait <- trait
  t1g$value <- "Tbv"
  t1g$Population <- "GN1"
  partitionDF1 <- rbind(partitionDF1, t1g)
  
  #Programme 2
  # EBV
  t2 <- Part2GNSummary[[paste0("Ebv", trait, "_s")]]$abs
  t2$Program <- "PN2"
  t2$Trait <- trait
  t2$value <- "Ebv"
  t2$Population <- "GN2"
  partitionDF2 <- rbind(partitionDF2, t2)
  
  # TBV
  t2g <- Part2gGNSummary[[paste0("Tbv", trait, "_s")]]$abs
  t2g$Program <- "PN2"
  t2g$Trait <- trait
  t2g$value <- "Tbv"
  t2g$Population <- "GN2"
  partitionDF2 <- rbind(partitionDF2, t2g)
  
  rm(t1, t1g, t2, t2g)
  
  # MULTIPLIER
  # Programme 1
  # EBV
  t1 <- Part1PNSummary[[paste0("Ebv", trait, "_s")]]$abs
  t1$Program <- "PN1"
  t1$Trait <- trait
  t1$value <- "Ebv"
  t1$Population <- "PN1"
  partitionDF1 <- rbind(partitionDF1, t1)
  
  # TBV
  t1g <- Part1gPNSummary[[paste0("Tbv", trait, "_s")]]$abs
  t1g$Program <- "PN1"
  t1g$Trait <- trait
  t1g$value <- "Tbv"
  t1g$Population <- "PN1"
  partitionDF1 <- rbind(partitionDF1, t1g)
  
  # Programme 2
  # EBV
  t2 <- Part2PNSummary[[paste0("Ebv", trait, "_s")]]$abs
  t2$Program <- "PN2"
  t2$Trait <- trait
  t2$value <- "Ebv"
  t2$Population <- "PN2"
  partitionDF2 <- rbind(partitionDF2, t2)
  
  #TBV
  t2g <- Part2gPNSummary[[paste0("Tbv", trait, "_s")]]$abs
  t2g$Program <- "PN2"
  t2g$Trait <- trait
  t2g$value <- "Tbv"
  t2g$Population <- "PN2"
  partitionDF2 <- rbind(partitionDF2, t2g)
  
}

#write partition results  
write.table(partitionDF1, paste0("PartitionPN1.csv"), quote=FALSE, row.names=FALSE)
write.table(partitionDF2, paste0("PartitionPN2.csv"), quote=FALSE, row.names=FALSE)



