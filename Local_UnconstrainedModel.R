###################################################################################
#
# Wisconsin bumble bee data
#
# Step 2: Urban/Agricultural Pools to Local Greenspaces
# Revised analysis with unconstrained null model (maintains species richness, but
# not species frequency)
#
# Community Trait Means
# Taxonomic Diversity
# Functional Diversity
# Phylogenetic Diversity
# Null Models
# PLS Analyses
#
# KI Perry; 27 September 2024
#
###################################################################################

## import species and trait data
t <- read.csv("./bb_rtraits_v3.csv", row.names=1)
a <- read.csv("./bb_abund_v3.csv", row.names=1)

str(a)
#change dataset to presence/absence
a[a > 0] <- 1
str(a)
rowSums(a)

# create a vector with the column sums for each species
# species not collected in Madison will have a 0
sp <- colSums(a)
sp

# removes any columns (i.e. species) that were note collected in Madison
a <- a[, colSums(a != 0) > 0]
str(a)

# add the sp vector as a column in the trait matrix, shows which species
# were collected in Madison and which were absent (i.e. with a 0)
t$sp <- sp
t

# use the sp values to remove rows of species not collected in Madison
# then remove the column because we don't need it anymore
t <- t[t$sp != 0, ]
t <- t[,-28]
str(t)

# trim the trait dataset by identifying traits that are highly correlated
# or lack sufficient variance among species
names(t)
str(t)
plot(t)
cor(t, method = c("pearson"), use = "complete.obs")

# lecty, nest construction, sociality, activity, parasitism, and pollen transport lack sufficient
# variance among bumble bee species - remove these traits
t2 <- t[,-8]#lecty
t2 <- t2[,-9]#nest construction
t2 <- t2[,-9]#sociality
t2 <- t2[,-9]#activity
t2 <- t2[,-9]#parasitism
t2 <- t2[,-9]#pollen transport

plot(t2)
cor(t2, method = c("pearson"), use = "complete.obs")

t2 <- t2[,-21]#wt

t2 <- t2[,-20]#corbicula width

t2 <- t2[,-10]#wing marginal cell length (keeping inter-tegular distance)
t2 <- t2[,-10]#wing width

t2 <- t2[,-12]#eye width

plot(t2)
cor(t2, method = c("pearson"), use = "complete.obs")

# remove the body size variables from the literature except body length variances
# these are correlated with each other and several other traits
t2 <- t2[,-1]#average queen body length
t2 <- t2[,-2]#average male body length
t2 <- t2[,-3]#average worker body length

plot(t2)
cor(t2, method = c("pearson"), use = "complete.obs")

t2 <- t2[,-10] #scape length
# will group other traits on head to limit their influence with functional diversity

#copy dataset in case transformataions are needed
t3 <- t2

str(t3)
t3$nestl <- as.factor(t3$nestl)
t3$tl <- as.factor(t3$tl)
str(t3)

#check traits for normality
hist(t2$qbl_var)
hist(log(t2$qbl_var))#
t3$qbl_var <- log(t3$qbl_var + 1)

hist(t2$mbl_var)
hist(log(t2$mbl_var))#
t3$mbl_var <- log(t3$mbl_var + 1)

hist(t2$wbl_var)
hist(log(t2$wbl_var))

hist(t2$it)
hist(log(t2$it))#
t3$it <- log(t3$it + 1)

hist(t2$rwingl)
hist(log(t2$rwingl))
t3$rwingl <- log(t3$rwingl + 1)

hist(t2$rheadw)
hist(log(t2$rheadw))

hist(t2$reyel)
hist(log(t2$reyel) + 1)#
t3$reyel <- log(t3$reyel + 1)

hist(t2$thairl)
hist(log(t2$thairl))#
t3$thairl <- log(t3$thairl + 1)

hist(t2$bsetael)
hist(log(t2$bsetael))#

hist(t2$rbtl)
hist(log(t2$rbtl))#
t3$rbtl <- log(t3$rbtl + 1)

# import phylogenetic tree
library(ape)
library(phytools)
library(phylotools)

p <- ape::read.nexus(file = "bbtree.nex")
plotTree(p, ftype="i")
Ntip(p)
p$tip.label

# trim tree to include only those species collected within Madison
names(a)
p.species <- c("boim","bobi","bova","boaf","bogr","boru","bofe","bobo","boau")
p2 <- drop.tip(p, setdiff(p$tip.label, p.species))
plotTree(p2, ftype="i")

# double check that all species are present in both datasets
intersect(colnames(a), rownames(t3))
intersect(colnames(a), p2$tip.label)
intersect(p2$tip.label, rownames(t3))
names(a)

# double check if a species is present in one dataset but not the other
setdiff(colnames(a), rownames(t3))
setdiff(rownames (t3), colnames(a))
setdiff(rownames (t3), p2$tip.label)
setdiff(colnames(a), p2$tip.label)

# double check all species names are in the same order
rownames(t3) == colnames(a) 
colnames(a) == p2$tip.label
rownames(t3) == p2$tip.label

# we are good to go!

##############################################################################
## Community Metrics

library(FD)
library(picante)
library(gawdis)
library(betapart)

####
#observed functional CWM
cwm.obs <- functcomp(t3, as.matrix(a), CWM.type = "all")
cwm.obs

####
# observed taxonomic beta diversity
# create beta part object for analyses
str(a)
bb.core <- betapart.core(a)

# returns three dissimilarity matrices containing 
# pairwise between-site values of each beta-diversity component
bb.dist <- beta.pair(bb.core, index.family = "sorensen")
str(bb.dist)

####
# observed functional beta diversity
# weight traits on the head together to limit their total influence on the metric
tdis <- gawdis(t3, w.type = "optimized", opti.maxiter = 300,
               groups.weight = T, groups = c(1, 2, 3, 4, 5, 6, 7, 8, 8, 9, 10, 11))
attr(tdis, "correls")
attr(tdis, "weights")

# calculate the distance matrix
pcoBB <- dudi.pco(sqrt(tdis), scannf = FALSE, nf = 4)#select four axes
scatter(pcoBB)

pcoBB$li
sum(pcoBB$eig[1:4]) / sum(pcoBB$eig)#0.74
sum(pcoBB$eig[1:3]) / sum(pcoBB$eig)#0.61
sum(pcoBB$eig[1:2]) / sum(pcoBB$eig)#0.43

# check correlations among axes and traits
str(t2)
cor(pcoBB$li, t2, use = "complete.obs")

# due to number of bee species at each site (one site has only 3 bee species),
# we can only use first two axes of PCoA for functional diversity metrics
t.ax <- as.matrix(pcoBB$li[1:2])

# returns pairwise between-site values of each functional beta-diversity component
bb.fun <- functional.beta.pair(a, t.ax, index.family = "sorensen")
str(bb.fun)

# observed functional alpha diversity
# run rao function first
bb.rao <- Rao(sample = t(a), dfunc = tdis, dphyl = NULL, weight = FALSE, Jost = TRUE, structure = NULL)
bb.falpha <- bb.rao$FD$Alpha
bb.falpha

####
# observed phylogenetic diversity
# create beta part object for analyses
bb.pcore <- phylo.betapart.core(a, p2)

# returns pairwise between-site values of each phylogenetic beta-diversity component
bb.phy <- phylo.beta.pair(bb.pcore, index.family = "sorensen")
str(bb.phy)

########################################################################################
# Null Model - Urban Pool to Local Greenspaces
# run the null model with 999 iterations
numberReps <- 999

# create empty matrices to store the results of each iteration of the null model:
# cwms
nqbl_var <- nmbl_var <- nwbl_var <- ntl_0 <- ntl_1 <- ntl_2 <- nnestl_0 <- nnestl_2 <- nit <- nwingl <- nheadw <- neyel <- nthairl <- nbsetael <- nbtl <- matrix(NA,
                                                                                                                                                                 nrow = nrow(a), ncol = numberReps, dimnames = list(rownames(a), paste0("n", 1:numberReps)))

# taxonomic beta diversity
nbsim <- nbsne <- nbsor <- matrix(NA, nrow = nrow(a), ncol = numberReps, 
                                  dimnames = list(rownames(a), paste0("n", 1:numberReps)))

# functional alpha and beta diversity
nfalpha <- nfsim <- nfsne <- nfsor <- matrix(NA, nrow = nrow(a), ncol = numberReps, 
                                             dimnames = list(rownames(a), paste0("n", 1:numberReps)))

# phylogenetic beta diversity
npsim <- npsne <- npsor <- matrix(NA, nrow = nrow(a), ncol = numberReps, 
                                  dimnames = list(rownames(a), paste0("n", 1:numberReps)))

#create null model for each repetition:

for(i in 1:numberReps){
  print(i) 
  
  # randomize trait matrix
  ntraits <- t3[sample(1:nrow(t3)),]
  rownames(ntraits) <- rownames(t3)
  
  # randomize presence/absence matrix
  # independent swap constrains by species richness and frequency
  #spBB <- randomizeMatrix(samp = a, null.model = "independentswap")
  spBB <- randomizeMatrix(samp = a, null.model = "richness") # Change in null model from independent swap to richness
  print(rownames(ntraits) == colnames(spBB)) 
  
  # randomize trait distance matrix
  ntdis <- gawdis(ntraits, w.type = "optimized", opti.maxiter = 300,
                  groups.weight = T, groups = c(1, 2, 3, 4, 5, 6, 7, 8, 8, 9, 10, 11))
  
  # randomize phylogenetic tree
  np <- tipShuffle(p2)
  
  # CWM calculations
  cwm.null <- functcomp(x = ntraits, a = as.matrix(spBB), CWM.type = "all")
  nqbl_var[,i] <- cwm.null$qbl_var
  nmbl_var[,i] <- cwm.null$mbl_var
  nwbl_var[,i] <- cwm.null$wbl_var
  ntl_0[,i] <- cwm.null$tl_0
  ntl_1[,i] <- cwm.null$tl_1
  ntl_2[,i] <- cwm.null$tl_2
  nnestl_0[,i] <- cwm.null$nestl_0
  nnestl_2[,i] <- cwm.null$nestl_2
  nit[,i] <- cwm.null$it
  nwingl[,i] <- cwm.null$rwingl
  nheadw[,i] <- cwm.null$rheadw
  neyel[,i] <- cwm.null$reyel
  nthairl[,i] <- cwm.null$thairl
  nbsetael[,i] <- cwm.null$bsetael
  nbtl[,i] <- cwm.null$rbtl
  
  # Functional alpha diversity
  nrao <- Rao(sample = t(spBB), dfunc = ntdis, dphyl = NULL, weight = FALSE, Jost = TRUE, structure = NULL)
  nfalpha[,i] <- nrao$FD$Alpha
  
  # Taxonomic beta diversity indices
  nbb.core <- betapart.core(spBB)
  nbb.dist <- beta.pair(nbb.core, index.family = "sorensen")
  nsim.dist <- as.matrix(nbb.dist$beta.sim)
  nsne.dist <- as.matrix(nbb.dist$beta.sne)
  nsor.dist <- as.matrix(nbb.dist$beta.sor)
  nbsim[,i] <- colMeans(nsim.dist)
  nbsne[,i] <- colMeans(nsne.dist)
  nbsor[,i] <- colMeans(nsor.dist)
  
  # Functional beta diversity indices
  npcoBB <- dudi.pco(sqrt(ntdis), scannf = FALSE, nf = 2)
  nt <- as.matrix(npcoBB$li)
  nbb.fun <- functional.beta.pair(spBB, nt, index.family = "sorensen")
  nfsim.dist <- as.matrix(nbb.fun$funct.beta.sim)
  nfsne.dist <- as.matrix(nbb.fun$funct.beta.sne)
  nfsor.dist <- as.matrix(nbb.fun$funct.beta.sor)
  nfsim[,i] <- colMeans(nfsim.dist)
  nfsne[,i] <- colMeans(nfsne.dist)
  nfsor[,i] <- colMeans(nfsor.dist)
  
  # Phylogenetic beta diversity indices
  mspBB <- match.phylo.comm(np, spBB)
  print(colnames(mspBB$comm) == np$tip.label)
  nbb.pcore <- phylo.betapart.core(mspBB$comm, np)
  nbb.phy <- phylo.beta.pair(nbb.pcore, index.family = "sorensen")
  npsim.dist <- as.matrix(nbb.phy$phylo.beta.sim)
  npsne.dist <- as.matrix(nbb.phy$phylo.beta.sne)
  npsor.dist <- as.matrix(nbb.phy$phylo.beta.sor)
  npsim[,i] <- colMeans(npsim.dist)
  npsne[,i] <- colMeans(npsne.dist)
  npsor[,i] <- colMeans(npsor.dist)
  
}

write.csv(nqbl_var, file = "Urban_to_Local_2Axes_Unconstrained/nqbl_var.u.csv")
write.csv(nmbl_var, file = "Urban_to_Local_2Axes_Unconstrained/nmbl_var.u.csv")
write.csv(nwbl_var, file = "Urban_to_Local_2Axes_Unconstrained/nwbl_var.u.csv")
write.csv(ntl_0, file = "Urban_to_Local_2Axes_Unconstrained/ntl_0.u.csv")
write.csv(ntl_1, file = "Urban_to_Local_2Axes_Unconstrained/ntl_1.u.csv")
write.csv(ntl_2, file = "Urban_to_Local_2Axes_Unconstrained/ntl_2.u.csv")
write.csv(nnestl_0, file = "Urban_to_Local_2Axes_Unconstrained/nnestl_0.u.csv")
write.csv(nnestl_2, file = "Urban_to_Local_2Axes_Unconstrained/nnestl_2.u.csv")
write.csv(nit, file = "Urban_to_Local_2Axes_Unconstrained/nit.u.csv")
write.csv(nwingl, file = "Urban_to_Local_2Axes_Unconstrained/nwingl.u.csv")
write.csv(nheadw, file = "Urban_to_Local_2Axes_Unconstrained/nheadw.u.csv")
write.csv(neyel, file = "Urban_to_Local_2Axes_Unconstrained/neyel.u.csv")
write.csv(nthairl, file = "Urban_to_Local_2Axes_Unconstrained/nthairl.u.csv")
write.csv(nbsetael, file = "Urban_to_Local_2Axes_Unconstrained/nbsetael.u.csv")
write.csv(nbtl, file = "Urban_to_Local_2Axes_Unconstrained/nbtl.u.csv")

write.csv(nbsim, file = "Urban_to_Local_2Axes_Unconstrained/bb_tbeta_sim.u.csv")
write.csv(nbsne, file = "Urban_to_Local_2Axes_Unconstrained/bb_tbeta_sne.u.csv")
write.csv(nbsor, file = "Urban_to_Local_2Axes_Unconstrained/bb_tbeta_sor.u.csv")

write.csv(nfalpha, file = "Urban_to_Local_2Axes_Unconstrained/bb_falpha.u.csv")
write.csv(nfsim, file = "Urban_to_Local_2Axes_Unconstrained/bb_fbeta_sim.u.csv")
write.csv(nfsne, file = "Urban_to_Local_2Axes_Unconstrained/bb_fbeta_sne.u.csv")
write.csv(nfsor, file = "Urban_to_Local_2Axes_Unconstrained/bb_fbeta_sor.u.csv")

write.csv(npsim, file = "Urban_to_Local_2Axes_Unconstrained/bb_pbeta_sim.csv")
write.csv(npsne, file = "Urban_to_Local_2Axes_Unconstrained/bb_pbeta_sne.csv")
write.csv(npsor, file = "Urban_to_Local_2Axes_Unconstrained/bb_pbeta_sor.csv")


########################################################################################################
# Now, let's run the model and not shuffle the trait matrix, only randomize the abundance matrix
# based on feedback from a reviewier

## import species and trait data
t <- read.csv("./bb_rtraits_v3.csv", row.names=1)
a <- read.csv("./bb_abund_v3.csv", row.names=1)

str(a)
#change dataset to presence/absence
a[a > 0] <- 1
str(a)
rowSums(a)

# create a vector with the column sums for each species
# species not collected in Madison will have a 0
sp <- colSums(a)
sp

# removes any columns (i.e. species) that were note collected in Madison
a <- a[, colSums(a != 0) > 0]
str(a)

# add the sp vector as a column in the trait matrix, shows which species
# were collected in Madison and which were absent (i.e. with a 0)
t$sp <- sp
t

# use the sp values to remove rows of species not collected in Madison
# then remove the column because we don't need it anymore
t <- t[t$sp != 0, ]
t <- t[,-28]
str(t)

# trim the trait dataset by identifying traits that are highly correlated
# or lack sufficient variance among species
names(t)
str(t)
plot(t)
cor(t, method = c("pearson"), use = "complete.obs")

# lecty, nest construction, sociality, activity, parasitism, and pollen transport lack sufficient
# variance among bumble bee species - remove these traits
t2 <- t[,-8]#lecty
t2 <- t2[,-9]#nest construction
t2 <- t2[,-9]#sociality
t2 <- t2[,-9]#activity
t2 <- t2[,-9]#parasitism
t2 <- t2[,-9]#pollen transport

plot(t2)
cor(t2, method = c("pearson"), use = "complete.obs")

t2 <- t2[,-21]#wt

t2 <- t2[,-20]#corbicula width

t2 <- t2[,-10]#wing marginal cell length (keeping inter-tegular distance)
t2 <- t2[,-10]#wing width

t2 <- t2[,-12]#eye width

plot(t2)
cor(t2, method = c("pearson"), use = "complete.obs")

# remove the body size variables from the literature except body length variances
# these are correlated with each other and several other traits
t2 <- t2[,-1]#average queen body length
t2 <- t2[,-2]#average male body length
t2 <- t2[,-3]#average worker body length

plot(t2)
cor(t2, method = c("pearson"), use = "complete.obs")

t2 <- t2[,-10] #scape length
# will group other traits on head to limit their influence with functional diversity

#copy dataset in case transformataions are needed
t3 <- t2

str(t3)
t3$nestl <- as.factor(t3$nestl)
t3$tl <- as.factor(t3$tl)
str(t3)

#check traits for normality
hist(t2$qbl_var)
hist(log(t2$qbl_var))#
t3$qbl_var <- log(t3$qbl_var + 1)

hist(t2$mbl_var)
hist(log(t2$mbl_var))#
t3$mbl_var <- log(t3$mbl_var + 1)

hist(t2$wbl_var)
hist(log(t2$wbl_var))

hist(t2$it)
hist(log(t2$it))#
t3$it <- log(t3$it + 1)

hist(t2$rwingl)
hist(log(t2$rwingl))
t3$rwingl <- log(t3$rwingl + 1)

hist(t2$rheadw)
hist(log(t2$rheadw))

hist(t2$reyel)
hist(log(t2$reyel) + 1)#
t3$reyel <- log(t3$reyel + 1)

hist(t2$thairl)
hist(log(t2$thairl))#
t3$thairl <- log(t3$thairl + 1)

hist(t2$bsetael)
hist(log(t2$bsetael))#

hist(t2$rbtl)
hist(log(t2$rbtl))#
t3$rbtl <- log(t3$rbtl + 1)

# import phylogenetic tree
library(ape)
library(phytools)
library(phylotools)

p <- ape::read.nexus(file = "bbtree.nex")
plotTree(p, ftype="i")
Ntip(p)
p$tip.label

# trim tree to include only those species collected within Madison
names(a)
p.species <- c("boim","bobi","bova","boaf","bogr","boru","bofe","bobo","boau")
p2 <- drop.tip(p, setdiff(p$tip.label, p.species))
plotTree(p2, ftype="i")

# double check that all species are present in both datasets
intersect(colnames(a), rownames(t3))
intersect(colnames(a), p2$tip.label)
intersect(p2$tip.label, rownames(t3))
names(a)

# double check if a species is present in one dataset but not the other
setdiff(colnames(a), rownames(t3))
setdiff(rownames (t3), colnames(a))
setdiff(rownames (t3), p2$tip.label)
setdiff(colnames(a), p2$tip.label)

# double check all species names are in the same order
rownames(t3) == colnames(a) 
colnames(a) == p2$tip.label
rownames(t3) == p2$tip.label

# we are good to go!

##############################################################################
## Community Metrics

library(FD)
library(picante)
library(gawdis)
library(betapart)

####
#observed functional CWM
cwm.obs <- functcomp(t3, as.matrix(a), CWM.type = "all")
cwm.obs

####
# observed taxonomic beta diversity
# create beta part object for analyses
str(a)
bb.core <- betapart.core(a)

# returns three dissimilarity matrices containing 
# pairwise between-site values of each beta-diversity component
bb.dist <- beta.pair(bb.core, index.family = "sorensen")
str(bb.dist)

####
# observed functional beta diversity
# weight traits on the head together to limit their total influence on the metric
tdis <- gawdis(t3, w.type = "optimized", opti.maxiter = 300,
               groups.weight = T, groups = c(1, 2, 3, 4, 5, 6, 7, 8, 8, 9, 10, 11))
attr(tdis, "correls")
attr(tdis, "weights")

# calculate the distance matrix
pcoBB <- dudi.pco(sqrt(tdis), scannf = FALSE, nf = 4)#select four axes
scatter(pcoBB)

pcoBB$li
sum(pcoBB$eig[1:4]) / sum(pcoBB$eig)#0.74
sum(pcoBB$eig[1:3]) / sum(pcoBB$eig)#0.61
sum(pcoBB$eig[1:2]) / sum(pcoBB$eig)#0.43

# check correlations among axes and traits
str(t2)
cor(pcoBB$li, t2, use = "complete.obs")

# due to number of bee species at each site (one site has only 3 bee species),
# we can only use first two axes of PCoA for functional diversity metrics
t.ax <- as.matrix(pcoBB$li[1:2])

# returns pairwise between-site values of each functional beta-diversity component
bb.fun <- functional.beta.pair(a, t.ax, index.family = "sorensen")
str(bb.fun)

# observed functional alpha diversity
# run rao function first
bb.rao <- Rao(sample = t(a), dfunc = tdis, dphyl = NULL, weight = FALSE, Jost = TRUE, structure = NULL)
bb.falpha <- bb.rao$FD$Alpha
bb.falpha

####
# observed phylogenetic diversity
# create beta part object for analyses
bb.pcore <- phylo.betapart.core(a, p2)

# returns pairwise between-site values of each phylogenetic beta-diversity component
bb.phy <- phylo.beta.pair(bb.pcore, index.family = "sorensen")
str(bb.phy)

########################################################################################
# Null Model - Urban Pool to Local Greenspaces
# run the null model with 999 iterations
numberReps <- 999

# create empty matrices to store the results of each iteration of the null model:
# cwms
nqbl_var <- nmbl_var <- nwbl_var <- ntl_0 <- ntl_1 <- ntl_2 <- nnestl_0 <- nnestl_2 <- nit <- nwingl <- nheadw <- neyel <- nthairl <- nbsetael <- nbtl <- matrix(NA,
                                                                                                                                                                 nrow = nrow(a), ncol = numberReps, dimnames = list(rownames(a), paste0("n", 1:numberReps)))

# taxonomic beta diversity
nbsim <- nbsne <- nbsor <- matrix(NA, nrow = nrow(a), ncol = numberReps, 
                                  dimnames = list(rownames(a), paste0("n", 1:numberReps)))

# functional alpha and beta diversity
nfalpha <- nfsim <- nfsne <- nfsor <- matrix(NA, nrow = nrow(a), ncol = numberReps, 
                                             dimnames = list(rownames(a), paste0("n", 1:numberReps)))

# phylogenetic beta diversity
npsim <- npsne <- npsor <- matrix(NA, nrow = nrow(a), ncol = numberReps, 
                                  dimnames = list(rownames(a), paste0("n", 1:numberReps)))

#create null model for each repetition:

for(i in 1:numberReps){
  print(i) 
  
  # randomize presence/absence matrix
  # independent swap constrains by species richness and frequency
  spBB <- randomizeMatrix(samp = a, null.model = "independentswap")
  print(rownames(t3) == colnames(spBB)) 
  
  # calculate trait distance matrix
  ntdis <- gawdis(t3, w.type = "optimized", opti.maxiter = 300,
                  groups.weight = T, groups = c(1, 2, 3, 4, 5, 6, 7, 8, 8, 9, 10, 11))
  
  # randomize phylogenetic tree
  np <- tipShuffle(p2)
  
  # CWM calculations
  cwm.null <- functcomp(x = t3, a = as.matrix(spBB), CWM.type = "all")
  nqbl_var[,i] <- cwm.null$qbl_var
  nmbl_var[,i] <- cwm.null$mbl_var
  nwbl_var[,i] <- cwm.null$wbl_var
  ntl_0[,i] <- cwm.null$tl_0
  ntl_1[,i] <- cwm.null$tl_1
  ntl_2[,i] <- cwm.null$tl_2
  nnestl_0[,i] <- cwm.null$nestl_0
  nnestl_2[,i] <- cwm.null$nestl_2
  nit[,i] <- cwm.null$it
  nwingl[,i] <- cwm.null$rwingl
  nheadw[,i] <- cwm.null$rheadw
  neyel[,i] <- cwm.null$reyel
  nthairl[,i] <- cwm.null$thairl
  nbsetael[,i] <- cwm.null$bsetael
  nbtl[,i] <- cwm.null$rbtl
  
  # Functional alpha diversity
  nrao <- Rao(sample = t(spBB), dfunc = ntdis, dphyl = NULL, weight = FALSE, Jost = TRUE, structure = NULL)
  nfalpha[,i] <- nrao$FD$Alpha
  
  # Taxonomic beta diversity indices
  nbb.core <- betapart.core(spBB)
  nbb.dist <- beta.pair(nbb.core, index.family = "sorensen")
  nsim.dist <- as.matrix(nbb.dist$beta.sim)
  nsne.dist <- as.matrix(nbb.dist$beta.sne)
  nsor.dist <- as.matrix(nbb.dist$beta.sor)
  nbsim[,i] <- colMeans(nsim.dist)
  nbsne[,i] <- colMeans(nsne.dist)
  nbsor[,i] <- colMeans(nsor.dist)
  
  # Functional beta diversity indices
  npcoBB <- dudi.pco(sqrt(ntdis), scannf = FALSE, nf = 2)
  nt <- as.matrix(npcoBB$li)
  nbb.fun <- functional.beta.pair(spBB, nt, index.family = "sorensen")
  nfsim.dist <- as.matrix(nbb.fun$funct.beta.sim)
  nfsne.dist <- as.matrix(nbb.fun$funct.beta.sne)
  nfsor.dist <- as.matrix(nbb.fun$funct.beta.sor)
  nfsim[,i] <- colMeans(nfsim.dist)
  nfsne[,i] <- colMeans(nfsne.dist)
  nfsor[,i] <- colMeans(nfsor.dist)
  
  # Phylogenetic beta diversity indices
  mspBB <- match.phylo.comm(np, spBB)
  print(colnames(mspBB$comm) == np$tip.label)
  nbb.pcore <- phylo.betapart.core(mspBB$comm, np)
  nbb.phy <- phylo.beta.pair(nbb.pcore, index.family = "sorensen")
  npsim.dist <- as.matrix(nbb.phy$phylo.beta.sim)
  npsne.dist <- as.matrix(nbb.phy$phylo.beta.sne)
  npsor.dist <- as.matrix(nbb.phy$phylo.beta.sor)
  npsim[,i] <- colMeans(npsim.dist)
  npsne[,i] <- colMeans(npsne.dist)
  npsor[,i] <- colMeans(npsor.dist)
  
}

write.csv(nqbl_var, file = "Urban_to_Local_2Axes_NoTraitShuffle/nqbl_var.u.csv")
write.csv(nmbl_var, file = "Urban_to_Local_2Axes_NoTraitShuffle/nmbl_var.u.csv")
write.csv(nwbl_var, file = "Urban_to_Local_2Axes_NoTraitShuffle/nwbl_var.u.csv")
write.csv(ntl_0, file = "Urban_to_Local_2Axes_NoTraitShuffle/ntl_0.u.csv")
write.csv(ntl_1, file = "Urban_to_Local_2Axes_NoTraitShuffle/ntl_1.u.csv")
write.csv(ntl_2, file = "Urban_to_Local_2Axes_NoTraitShuffle/ntl_2.u.csv")
write.csv(nnestl_0, file = "Urban_to_Local_2Axes_NoTraitShuffle/nnestl_0.u.csv")
write.csv(nnestl_2, file = "Urban_to_Local_2Axes_NoTraitShuffle/nnestl_2.u.csv")
write.csv(nit, file = "Urban_to_Local_2Axes_NoTraitShuffle/nit.u.csv")
write.csv(nwingl, file = "Urban_to_Local_2Axes_NoTraitShuffle/nwingl.u.csv")
write.csv(nheadw, file = "Urban_to_Local_2Axes_NoTraitShuffle/nheadw.u.csv")
write.csv(neyel, file = "Urban_to_Local_2Axes_NoTraitShuffle/neyel.u.csv")
write.csv(nthairl, file = "Urban_to_Local_2Axes_NoTraitShuffle/nthairl.u.csv")
write.csv(nbsetael, file = "Urban_to_Local_2Axes_NoTraitShuffle/nbsetael.u.csv")
write.csv(nbtl, file = "Urban_to_Local_2Axes_NoTraitShuffle/nbtl.u.csv")

write.csv(nbsim, file = "Urban_to_Local_2Axes_NoTraitShuffle/bb_tbeta_sim.u.csv")
write.csv(nbsne, file = "Urban_to_Local_2Axes_NoTraitShuffle/bb_tbeta_sne.u.csv")
write.csv(nbsor, file = "Urban_to_Local_2Axes_NoTraitShuffle/bb_tbeta_sor.u.csv")

write.csv(nfalpha, file = "Urban_to_Local_2Axes_NoTraitShuffle/bb_falpha.u.csv")
write.csv(nfsim, file = "Urban_to_Local_2Axes_NoTraitShuffle/bb_fbeta_sim.u.csv")
write.csv(nfsne, file = "Urban_to_Local_2Axes_NoTraitShuffle/bb_fbeta_sne.u.csv")
write.csv(nfsor, file = "Urban_to_Local_2Axes_NoTraitShuffle/bb_fbeta_sor.u.csv")

write.csv(npsim, file = "Urban_to_Local_2Axes_NoTraitShuffle/bb_pbeta_sim.csv")
write.csv(npsne, file = "Urban_to_Local_2Axes_NoTraitShuffle/bb_pbeta_sne.csv")
write.csv(npsor, file = "Urban_to_Local_2Axes_NoTraitShuffle/bb_pbeta_sor.csv")
