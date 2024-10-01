###################################################################################
#
# Wisconsin bumble bee data
#
# Step 1: Regional Pool to Urban/Agricultural Pools
# Revised analysis that removed the trait shuffle randomization from the null model
#
# Community Trait Means
# Taxonomic Diversity
# Functional Diversity
# Phylogenetic Diversity
# Null Models
#
# KI Perry; 30 September 2024
#
###################################################################################

## Species data
a <- read.csv("./bb_abund_v3.csv", row.names=1)
str(a)
#change dataset to presence/absence for this part of the analyses
a[a > 0] <- 1
str(a)
rowSums(a)
colSums(a)


## Trait data
t <- read.csv("./bb_rtraits_v3.csv", row.names=1)
#trim the trait dataset by identifying traits that are highly correlated
#or lack sufficient variance among species
names(t)
str(t)
plot(t)

#lecty, nest construction, sociality, activity, parasitism, and pollen transport lack sufficient
#variance among bumble bee species - remove these traits
t2 <- t[,-8]#lecty
t2 <- t2[,-9]#nest construction
t2 <- t2[,-9]#sociality
t2 <- t2[,-9]#activity
t2 <- t2[,-9]#parasitism
t2 <- t2[,-9]#pollen transport

plot(t2)
cor(t2, method = c("pearson"), use = "complete.obs")

#remove the body size variables from the literature except body length variances
#these are correlated with each other and several other traits
t2 <- t2[,-1]#average queen body length
t2 <- t2[,-2]#average male body length
t2 <- t2[,-3]#average worker body length

#highly correlated with other traits
t2 <- t2[,-18]#wt

plot(t2)

library(ggplot2)
library(GGally)

cp <- ggpairs(t2, upper = list(continuous = wrap("cor", size = 5, color = "black")))
cp + theme(strip.text.x = element_text(size = 18), strip.text.y = element_text(size = 10))
cor(t2, method = c("pearson"), use = "complete.obs")

# Pairplot - for Appendix
png("Corplot_Regional.png", width = 1800, height = 1000, pointsize = 40)
cp <- ggpairs(t2, upper = list(continuous = wrap("cor", size = 5, color = "black")))
cp + theme(strip.text.x = element_text(size = 23), strip.text.y = element_text(size = 12))
dev.off()

t2 <- t2[,-17]#corbicula width

t2 <- t2[,-7]#wing marginal cell length (keeping inter-tegular distance)
t2 <- t2[,-7]#wing width

t2 <- t2[,-9]#eye width

t2 <- t2[,-10] #scape length
#will group other traits on head to limit their influence with functional diversity

ggpairs(t2)

#copy dataset in case transformations are needed
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
t3$it <- log(t3$it)

hist(t2$rwingl)
hist(log(t2$rwingl))

hist(t2$rheadw)
hist(log(t2$rheadw))

hist(t2$reyel)
hist(log(t2$reyel) + 1)#
t3$reyel <- log(t3$reyel + 1)

hist(t2$thairl)
hist(log(t2$thairl))#
t3$thairl <- log(t3$thairl)

hist(t2$bsetael)
hist(log(t2$bsetael))#
t3$bsetael <- log(t3$bsetael + 1)

hist(t2$rbtl)
hist(log(t2$rbtl))#
t3$rbtl <- log(t3$rbtl + 1)

## Phylogenetic data
library(ape)
citation("ape")
library(phytools)
library(phylotools)

#import bumble bee phylogenetic tree (nexus file)
#skip this and import trimmed tree below
#p <- ape::read.nexus(file = "combdivBayesiv19v3.nex.con")
#class(p)
#str(p)

#plotTree(p, ftype = "i", fsize = 0.5, lwd = 1)
#Ntip(p)

#check species in tree
#p$con_50_majrule$tip.label
#p.species <- c("perplexus","ternarius","impatiens","bimaculat",
#               "vagans","affinis","griseocol","rufocinct","fervidus","pensylv2",
#               "citrinus","borealis","auricomus")


#p2 <- as.phylo(p$con_50_majrule)
#p3 <- drop.tip(p2, setdiff(p2$tip.label, p.species))
#plotTree(p3, ftype="i")

#replace taxa labels so they match with those in the trait and presence/absence data
#n.species <- c("bope","botr","boim","bobi","bova","boaf","bogr","boru","bofe",
#               "bopn","boci","bobo","boau")

#bnames <- as.data.frame(cbind(p.species,n.species))

#p4 <- sub.taxa.label(p3, bnames)
#plotTree(p4, ftype="i")

#save revised tree
#ape::write.nexus(p4, file='bbtree.nex')

#import revised tree
p4 <- ape::read.nexus(file = "bbtree.nex")
plotTree(p4, ftype="i")
Ntip(p4)
p4$tip.label

#Double check that all species are present in both datasets
intersect(colnames(a), rownames(t3))
intersect(colnames(a), p4$tip.label)
intersect(p4$tip.label, rownames(t3))
names(a)

#Double check if a species is present in one dataset but not the other
setdiff(colnames(a), rownames(t3))
setdiff(rownames (t3), colnames(a))
setdiff(rownames (t3), p4$tip.label)
setdiff(colnames(a), p4$tip.label)

#Double check all species names are in the same order
rownames(t3) == colnames(a) 
colnames(a) == p4$tip.label
rownames(t3) == p4$tip.label

#We are good to go!

##############################################################################
## Observed Community Metrics

library(FD)
library(picante)
library(gawdis)
library(betapart)
library(dplyr)
citation(package = "ade4")
citation(package = "phytools")
citation(package = "picante")

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
# calculate the distance matrix
str(t3)
tdis <- gawdis(t3, w.type = "optimized", opti.maxiter = 300,
               groups.weight = T, groups = c(1, 2, 3, 4, 5, 6, 7, 8, 8, 9, 10, 11))
attr(tdis, "correls")
attr(tdis, "weights")

pcoBB <- dudi.pco(sqrt(tdis), scannf = FALSE, nf = 4)#select four axes
scatter(pcoBB)

pcoBB$li
sum(pcoBB$eig[1:4]) / sum(pcoBB$eig)#0.61
sum(pcoBB$eig[1:3]) / sum(pcoBB$eig)#0.50
sum(pcoBB$eig[1:2]) / sum(pcoBB$eig)#0.38

# check correlations among axes and traits
str(t2)
cor(pcoBB$li, t2, use = "complete.obs")

# due to number of bee species at each site (one site has only 3 bee species),
# we can only use first two axes of PCoA for functional diversity metrics
t.ax <- as.matrix(pcoBB$li[1:2])

# returns pairwise between-site values of each functional beta-diversity component
bb.fun <- functional.beta.pair(a, t.ax, index.family = "sorensen")
str(bb.fun)

####
# observed phylogenetic diversity
# create beta part object for analyses
bb.pcore <- phylo.betapart.core(a, p4)

# returns pairwise between-site values of each phylogenetic beta-diversity component
bb.phy <- phylo.beta.pair(bb.pcore, index.family = "sorensen")
str(bb.phy)

# calculate phylogenetic signals for traits
# Blomberg's K (Blomberg, Garland, and Ives 2003)

# create an empty dataframe to store the values for each trait
physig <- matrix(NA, nrow = 1, ncol = 10)
colnames(physig) <- c("qbl_var", "mbl_var", "wbl_var", "it", "rwingl", "rheadw", "reyel", "thairl",
                      "bsetael", "rbtl")

physig[,1] <- phylosig(p4, t3$qbl_var, method = "K")
physig[,2] <- phylosig(p4, t3$mbl_var, method = "K")
physig[,3] <- phylosig(p4, t3$wbl_var, method = "K")
physig[,4] <- phylosig(p4, t3$it, method = "K")
physig[,5] <- phylosig(p4, t3$rwingl, method = "K")
physig[,6] <- phylosig(p4, t3$rheadw, method = "K")
physig[,7] <- phylosig(p4, t3$reyel, method = "K")
physig[,8] <- phylosig(p4, t3$thairl, method = "K")
physig[,9] <- phylosig(p4, t3$bsetael, method = "K")
physig[,10] <- phylosig(p4, t3$rbtl, method = "K")

# needs to be a data frame for use later
physig <- as.data.frame(physig)

################################################################################
# test randomization code to make sure it is working correctly
# single iteration of the null model

# randomize community data matrix within samples (maintains sample species richness but not species frequency)
# because only fixing species richness, susceptible to type 1 errors
# some species identified in the regional pool were not collected, so for this model, we cannot
# fix patterns of species frequency

#for (i in 1:nrow(a)){
#  randomizedBB <- randomizeMatrix(samp = a, null.model = "richness", iterations = 1)
#}

#colSums(a)
#colSums(randomizedBB) #it worked!

#null CWMs
#cwm.test <- functcomp(t3, as.matrix(randomizedBB), CWM.type = "all")
#cwm.test

# randomize the trait matrix by randomizing rows to maintain trait interrelationships
# just changes the species name associated with a row of trait values

#traitsRand <- t3[sample(1:nrow(t3)),]
#rownames(traitsRand) <- rownames(t3)

# randomize the species names on the phylogenetic tree
#phyrand <- tipShuffle(p4)
#plotTree(phyrand, ftype="i")

# check code to match order of species in abundance and trait datasets with the randomized phylogeny
# diversity metrics cannot be calculated unless species are in the same order
#ma <- match.phylo.comm(phyrand, a)
#ma$comm
#mt <- match.phylo.data(phyrand, t3)
#mnt <- as.data.frame(mt$data)
#mnt <- within(mnt, {
#  qbl_var <- as.numeric(qbl_var)
#})
#str(mnt)
#phylosig(phyrand, mnt$qbl_var, method = "K")
#rownames(mnt) == phyrand$tip.label


########################################################################################
# Null Model - Regional Pool to Urban Pool
# run the null model with 999 iterations
numberReps <- 999

#create empty matrices to store the results of each iteration of the null model:
# cwms
nqbl_var <- nmbl_var <- nwbl_var <- ntl_0 <- ntl_1 <- ntl_2 <- nnestl_0 <- nnestl_1 <- nnestl_2 <- nit <- nwingl <- nheadw <- neyel <- nthairl <- nbsetael <- nbtl <- matrix(NA,
                                                                                                                                                                             nrow = nrow(a), ncol = numberReps, dimnames = list(rownames(a), paste0("n", 1:numberReps)))

# taxonomic beta diversity
nbsim <- nbsne <- nbsor <- matrix(NA, nrow = nrow(a), ncol = numberReps, 
                                  dimnames = list(rownames(a), paste0("n", 1:numberReps)))

# functional beta diversity
nfsim <- nfsne <- nfsor <- matrix(NA, nrow = nrow(a), ncol = numberReps, 
                                  dimnames = list(rownames(a), paste0("n", 1:numberReps)))

# phylogenetic beta diversity
npsim <- npsne <- npsor <- matrix(NA, nrow = nrow(a), ncol = numberReps, 
                                  dimnames = list(rownames(a), paste0("n", 1:numberReps)))

# phylogenetic signal
npsig <- matrix(NA, nrow = numberReps, ncol = ncol(t3), dimnames = list(paste0("n", 1:numberReps), colnames(t3)))

# create null model for each repetition:

for(i in 1:numberReps){
  print(i) 
  
  # randomize presence/absence matrix
  # richness constrains by species richness only
  spBB <- randomizeMatrix(samp = a, null.model = "richness")
  print(rownames(t3) == colnames(spBB)) 
  
  # randomize trait distance matrix
  ntdis <- gawdis(t3, w.type = "optimized", opti.maxiter = 300,
                  groups.weight = T, groups = c(1, 2, 3, 4, 5, 6, 7, 8, 8, 9, 10, 11))
  
  # randomize phylogenetic tree
  np <- tipShuffle(p4)
  
  # CWM calculations
  cwm.null <- functcomp(x = t3, a = as.matrix(spBB), CWM.type = "all")
  nqbl_var[,i] <- cwm.null$qbl_var
  nmbl_var[,i] <- cwm.null$mbl_var
  nwbl_var[,i] <- cwm.null$wbl_var
  ntl_0[,i] <- cwm.null$tl_0
  ntl_1[,i] <- cwm.null$tl_1
  ntl_2[,i] <- cwm.null$tl_2
  nnestl_0[,i] <- cwm.null$nestl_0
  nnestl_1[,i] <- cwm.null$nestl_1
  nnestl_2[,i] <- cwm.null$nestl_2
  nit[,i] <- cwm.null$it
  nwingl[,i] <- cwm.null$rwingl
  nheadw[,i] <- cwm.null$rheadw
  neyel[,i] <- cwm.null$reyel
  nthairl[,i] <- cwm.null$thairl
  nbsetael[,i] <- cwm.null$bsetael
  nbtl[,i] <- cwm.null$rbtl
  
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

# save the output files
write.csv(nqbl_var, file = "Regional_to_Urban_2Axes_NoTraitShuffle/nqbl_var.csv")
write.csv(nmbl_var, file = "Regional_to_Urban_2Axes_NoTraitShuffle/nmbl_var.csv")
write.csv(nwbl_var, file = "Regional_to_Urban_2Axes_NoTraitShuffle/nwbl_var.csv")
write.csv(ntl_0, file = "Regional_to_Urban_2Axes_NoTraitShuffle/ntl_0.csv")
write.csv(ntl_1, file = "Regional_to_Urban_2Axes_NoTraitShuffle/ntl_1.csv")
write.csv(ntl_2, file = "Regional_to_Urban_2Axes_NoTraitShuffle/ntl_2.csv")
write.csv(nnestl_0, file = "Regional_to_Urban_2Axes_NoTraitShuffle/nnestl_0.csv")
write.csv(nnestl_1, file = "Regional_to_Urban_2Axes_NoTraitShuffle/nnestl_1.csv")
write.csv(nnestl_2, file = "Regional_to_Urban_2Axes_NoTraitShuffle/nnestl_2.csv")
write.csv(nit, file = "Regional_to_Urban_2Axes_NoTraitShuffle/nit.csv")
write.csv(nwingl, file = "Regional_to_Urban_2Axes_NoTraitShuffle/nwingl.csv")
write.csv(nheadw, file = "Regional_to_Urban_2Axes_NoTraitShuffle/nheadw.csv")
write.csv(neyel, file = "Regional_to_Urban_2Axes_NoTraitShuffle/neyel.csv")
write.csv(nthairl, file = "Regional_to_Urban_2Axes_NoTraitShuffle/nthairl.csv")
write.csv(nbsetael, file = "Regional_to_Urban_2Axes_NoTraitShuffle/nbsetael.csv")
write.csv(nbtl, file = "Regional_to_Urban_2Axes_NoTraitShuffle/nbtl.csv")

write.csv(nbsim, file = "Regional_to_Urban_2Axes_NoTraitShuffle/bb_tbeta_sim.csv")
write.csv(nbsne, file = "Regional_to_Urban_2Axes_NoTraitShuffle/bb_tbeta_sne.csv")
write.csv(nbsor, file = "Regional_to_Urban_2Axes_NoTraitShuffle/bb_tbeta_sor.csv")

write.csv(nfsim, file = "Regional_to_Urban_2Axes_NoTraitShuffle/bb_fbeta_sim.csv")
write.csv(nfsne, file = "Regional_to_Urban_2Axes_NoTraitShuffle/bb_fbeta_sne.csv")
write.csv(nfsor, file = "Regional_to_Urban_2Axes_NoTraitShuffle/bb_fbeta_sor.csv")

write.csv(npsim, file = "Regional_to_Urban_2Axes_NoTraitShuffle/bb_pbeta_sim.csv")
write.csv(npsne, file = "Regional_to_Urban_2Axes_NoTraitShuffle/bb_pbeta_sne.csv")
write.csv(npsor, file = "Regional_to_Urban_2Axes_NoTraitShuffle/bb_pbeta_sor.csv")


# load the datasets
nqbl_var <- read.csv("Regional_to_Urban_2Axes_NoTraitShuffle/nqbl_var.csv", row.names=1)
nmbl_var <- read.csv("Regional_to_Urban_2Axes_NoTraitShuffle/nmbl_var.csv", row.names=1)
nwbl_var <- read.csv("Regional_to_Urban_2Axes_NoTraitShuffle/nwbl_var.csv", row.names=1)
ntl_0 <- read.csv("Regional_to_Urban_2Axes_NoTraitShuffle/ntl_0.csv", row.names=1)
ntl_1 <- read.csv("Regional_to_Urban_2Axes_NoTraitShuffle/ntl_1.csv", row.names=1)
ntl_2 <- read.csv("Regional_to_Urban_2Axes_NoTraitShuffle/ntl_2.csv", row.names=1)
nnestl_0 <- read.csv("Regional_to_Urban_2Axes_NoTraitShuffle/nnestl_0.csv", row.names=1)
nnestl_1 <- read.csv("Regional_to_Urban_2Axes_NoTraitShuffle/nnestl_1.csv", row.names=1)
nnestl_2 <- read.csv("Regional_to_Urban_2Axes_NoTraitShuffle/nnestl_2.csv", row.names=1)
nit <- read.csv("Regional_to_Urban_2Axes_NoTraitShuffle/nit.csv", row.names=1)
nwingl <- read.csv("Regional_to_Urban_2Axes_NoTraitShuffle/nwingl.csv", row.names=1)
nheadw <- read.csv("Regional_to_Urban_2Axes_NoTraitShuffle/nheadw.csv", row.names=1)
neyel <- read.csv("Regional_to_Urban_2Axes_NoTraitShuffle/neyel.csv", row.names=1)
nthairl <- read.csv("Regional_to_Urban_2Axes_NoTraitShuffle/nthairl.csv", row.names=1)
nbsetael <- read.csv("Regional_to_Urban_2Axes_NoTraitShuffle/nbsetael.csv", row.names=1)
nbtl <- read.csv("Regional_to_Urban_2Axes_NoTraitShuffle/nbtl.csv", row.names=1)

nbsim <- read.csv("Regional_to_Urban_2Axes_NoTraitShuffle/bb_tbeta_sim.csv", row.names=1)
nbsne <- read.csv("Regional_to_Urban_2Axes_NoTraitShuffle/bb_tbeta_sne.csv", row.names=1)
nbsor <- read.csv("Regional_to_Urban_2Axes_NoTraitShuffle/bb_tbeta_sor.csv", row.names=1)

nfsim <- read.csv("Regional_to_Urban_2Axes_NoTraitShuffle/bb_fbeta_sim.csv", row.names=1)
nfsne <- read.csv("Regional_to_Urban_2Axes_NoTraitShuffle/bb_fbeta_sne.csv", row.names=1)
nfsor <- read.csv("Regional_to_Urban_2Axes_NoTraitShuffle/bb_fbeta_sor.csv", row.names=1)

npsim <- read.csv("Regional_to_Urban_2Axes_NoTraitShuffle/bb_pbeta_sim.csv", row.names=1)
npsne <- read.csv("Regional_to_Urban_2Axes_NoTraitShuffle/bb_pbeta_sne.csv", row.names=1)
npsor <- read.csv("Regional_to_Urban_2Axes_NoTraitShuffle/bb_pbeta_sor.csv", row.names=1)
npsig <- read.csv("Regional_to_Urban_2Axes_NoTraitShuffle/bb_phylo_signal.csv", row.names=1)