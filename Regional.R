###################################################################################
#
# Wisconsin bumble bee data
#
# Step 1: Regional Pool to Urban/Agricultural Pools
#
# Community Trait Means
# Taxonomic Diversity
# Functional Diversity
# Phylogenetic Diversity
# Null Models
#
# KI Perry; 20 December 2021
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
  
  # randomize trait matrix
  ntraits <- t3[sample(1:nrow(t3)),]
  rownames(ntraits) <- rownames(t3)
  
  # randomize presence/absence matrix
  # richness constrains by species richness only
  spBB <- randomizeMatrix(samp = a, null.model = "richness")
  print(rownames(ntraits) == colnames(spBB)) 
  
  # randomize trait distance matrix
  ntdis <- gawdis(ntraits, w.type = "optimized", opti.maxiter = 300,
                  groups.weight = T, groups = c(1, 2, 3, 4, 5, 6, 7, 8, 8, 9, 10, 11))
  
  # randomize phylogenetic tree
  np <- tipShuffle(p4)
  
  # CWM calculations
  cwm.null <- functcomp(x = ntraits, a = as.matrix(spBB), CWM.type = "all")
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
  
  # Phylogenetic signal: Bloomberg's K
  mntraits <- match.phylo.data(np, ntraits)
  mnt <- data.frame(mntraits$data)
  mnt <- mnt %>% mutate_if(is.character,as.numeric)
  print(rownames(mnt) == np$tip.label)
  npsig[i,1] <- phylosig(np, mnt$qbl_var, method = "K")
  npsig[i,2] <- phylosig(np, mnt$mbl_var, method = "K")
  npsig[i,3] <- phylosig(np, mnt$wbl_var, method = "K")
  npsig[i,6] <- phylosig(np, mnt$it, method = "K")
  npsig[i,7] <- phylosig(np, mnt$rwingl, method = "K")
  npsig[i,8] <- phylosig(np, mnt$rheadw, method = "K")
  npsig[i,9] <- phylosig(np, mnt$reyel, method = "K")
  npsig[i,10] <- phylosig(np, mnt$thairl, method = "K")
  npsig[i,11] <- phylosig(np, mnt$bsetael, method = "K")
  npsig[i,12] <- phylosig(np, mnt$rbtl, method = "K")
  
}



write.csv(nqbl_var, file = "Regional_to_Urban_2Axes/nqbl_var.csv")
write.csv(nmbl_var, file = "Regional_to_Urban_2Axes/nmbl_var.csv")
write.csv(nwbl_var, file = "Regional_to_Urban_2Axes/nwbl_var.csv")
write.csv(ntl_0, file = "Regional_to_Urban_2Axes/ntl_0.csv")
write.csv(ntl_1, file = "Regional_to_Urban_2Axes/ntl_1.csv")
write.csv(ntl_2, file = "Regional_to_Urban_2Axes/ntl_2.csv")
write.csv(nnestl_0, file = "Regional_to_Urban_2Axes/nnestl_0.csv")
write.csv(nnestl_1, file = "Regional_to_Urban_2Axes/nnestl_1.csv")
write.csv(nnestl_2, file = "Regional_to_Urban_2Axes/nnestl_2.csv")
write.csv(nit, file = "Regional_to_Urban_2Axes/nit.csv")
write.csv(nwingl, file = "Regional_to_Urban_2Axes/nwingl.csv")
write.csv(nheadw, file = "Regional_to_Urban_2Axes/nheadw.csv")
write.csv(neyel, file = "Regional_to_Urban_2Axes/neyel.csv")
write.csv(nthairl, file = "Regional_to_Urban_2Axes/nthairl.csv")
write.csv(nbsetael, file = "Regional_to_Urban_2Axes/nbsetael.csv")
write.csv(nbtl, file = "Regional_to_Urban_2Axes/nbtl.csv")

write.csv(nbsim, file = "Regional_to_Urban_2Axes/bb_tbeta_sim.csv")
write.csv(nbsne, file = "Regional_to_Urban_2Axes/bb_tbeta_sne.csv")
write.csv(nbsor, file = "Regional_to_Urban_2Axes/bb_tbeta_sor.csv")

write.csv(nfsim, file = "Regional_to_Urban_2Axes/bb_fbeta_sim.csv")
write.csv(nfsne, file = "Regional_to_Urban_2Axes/bb_fbeta_sne.csv")
write.csv(nfsor, file = "Regional_to_Urban_2Axes/bb_fbeta_sor.csv")

write.csv(npsim, file = "Regional_to_Urban_2Axes/bb_pbeta_sim.csv")
write.csv(npsne, file = "Regional_to_Urban_2Axes/bb_pbeta_sne.csv")
write.csv(npsor, file = "Regional_to_Urban_2Axes/bb_pbeta_sor.csv")
write.csv(npsig, file = "Regional_to_Urban_2Axes/bb_phylo_signal.csv")


# load the datasets
nqbl_var <- read.csv("Regional_to_Urban_2Axes/nqbl_var.csv", row.names=1)
nmbl_var <- read.csv("Regional_to_Urban_2Axes/nmbl_var.csv", row.names=1)
nwbl_var <- read.csv("Regional_to_Urban_2Axes/nwbl_var.csv", row.names=1)
ntl_0 <- read.csv("Regional_to_Urban_2Axes/ntl_0.csv", row.names=1)
ntl_1 <- read.csv("Regional_to_Urban_2Axes/ntl_1.csv", row.names=1)
ntl_2 <- read.csv("Regional_to_Urban_2Axes/ntl_2.csv", row.names=1)
nnestl_0 <- read.csv("Regional_to_Urban_2Axes/nnestl_0.csv", row.names=1)
nnestl_1 <- read.csv("Regional_to_Urban_2Axes/nnestl_1.csv", row.names=1)
nnestl_2 <- read.csv("Regional_to_Urban_2Axes/nnestl_2.csv", row.names=1)
nit <- read.csv("Regional_to_Urban_2Axes/nit.csv", row.names=1)
nwingl <- read.csv("Regional_to_Urban_2Axes/nwingl.csv", row.names=1)
nheadw <- read.csv("Regional_to_Urban_2Axes/nheadw.csv", row.names=1)
neyel <- read.csv("Regional_to_Urban_2Axes/neyel.csv", row.names=1)
nthairl <- read.csv("Regional_to_Urban_2Axes/nthairl.csv", row.names=1)
nbsetael <- read.csv("Regional_to_Urban_2Axes/nbsetael.csv", row.names=1)
nbtl <- read.csv("Regional_to_Urban_2Axes/nbtl.csv", row.names=1)

nbsim <- read.csv("Regional_to_Urban_2Axes/bb_tbeta_sim.csv", row.names=1)
nbsne <- read.csv("Regional_to_Urban_2Axes/bb_tbeta_sne.csv", row.names=1)
nbsor <- read.csv("Regional_to_Urban_2Axes/bb_tbeta_sor.csv", row.names=1)

nfsim <- read.csv("Regional_to_Urban_2Axes/bb_fbeta_sim.csv", row.names=1)
nfsne <- read.csv("Regional_to_Urban_2Axes/bb_fbeta_sne.csv", row.names=1)
nfsor <- read.csv("Regional_to_Urban_2Axes/bb_fbeta_sor.csv", row.names=1)

npsim <- read.csv("Regional_to_Urban_2Axes/bb_pbeta_sim.csv", row.names=1)
npsne <- read.csv("Regional_to_Urban_2Axes/bb_pbeta_sne.csv", row.names=1)
npsor <- read.csv("Regional_to_Urban_2Axes/bb_pbeta_sor.csv", row.names=1)
npsig <- read.csv("Regional_to_Urban_2Axes/bb_phylo_signal.csv", row.names=1)

## Pval and SES:
#calculate standardized effect sizes (SES) for each trait and diversity metric
#the effect size is the difference between the observed value and the expected/ null values
#then divide the effect size by the standard deviation of the null distribution to get the standardized effect size
#this metric allows comparison among sites with different numbers of species

#test different methods with intertegular distance
#tested and they all work
#1
#meanNull_it <- rowMeans(nit)
#ES_it <- cwm.obs$it - meanNull_it
#sdNull_it <- apply(nit, 1, sd)
#SES_it.1 <- ES_it / sdNull_it

#2
#SES_it.2 <-(cwm.obs$it - rowMeans(nit)) / apply(nit, 1, sd, na.rm=T)

#3
#SES_it.3 <-(cwm.obs$it - apply(nit, MARGIN = 1, mean)) / apply(nit, MARGIN = 1, sd, na.rm=T)

#data.frame(SES_it.1, SES_it.2, SES_it.3) #all three produced the same results!

#################################################################################################################

# Community weighted means
cwm.obs.m <- as.matrix(t(colMeans(cwm.obs)))

## queen body length variance
SES_qblv <- (cwm.obs$qbl_var - apply(nqbl_var, MARGIN = 1, mean)) / apply(nqbl_var, MARGIN = 1, sd, na.rm=T)
SES_qblv

boxplot(SES_qblv)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

plot(SES_qblv, pch = 19, cex = 1.5)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

hist(as.matrix(nqbl_var), main = "Distribution of expected values - Queen Body Length Variance")
abline(v = cwm.obs.m[,1], col = "blue", lwd = 2)

pval.qblv <- apply(cbind(cwm.obs$qbl_var, nqbl_var), MARGIN = 1, rank)[1,] / 1000
pval.qblv
qblv <- cbind(cwm.obs$qbl_var, nqbl_var)
qblv.m <- as.matrix(t(colMeans(qblv)))
pval.qblv.a <- apply(qblv.m, MARGIN = 1, rank)[1,] / 1000
pval.qblv.a

w.qbl_var <- wilcox.test(SES_qblv, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)#
w.qbl_var


## male body length variance
SES_mblv <- (cwm.obs$mbl_var - apply(nmbl_var, MARGIN = 1, mean)) / apply(nmbl_var, MARGIN = 1, sd, na.rm=T)
SES_mblv

boxplot(SES_mblv)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

plot(SES_mblv, pch = 19, cex = 1.5)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

hist(as.matrix(nmbl_var), main = "Distribution of expected values - Male Body Length Variance")
abline(v = cwm.obs.m[,2], col = "blue", lwd = 2)

pval.mblv <- apply(cbind(cwm.obs$mbl_var, nmbl_var), MARGIN = 1, rank)[1,] / 1000
pval.mblv
mblv <- cbind(cwm.obs$mbl_var, nmbl_var)
mblv.m <- as.matrix(t(colMeans(mblv)))
pval.mblv.a <- apply(mblv.m, MARGIN = 1, rank)[1,] / 1000
pval.mblv.a

w.mbl_var <- wilcox.test(SES_mblv, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)#
w.mbl_var


## worker body length variance
SES_wblv <- (cwm.obs$wbl_var - apply(nwbl_var, MARGIN = 1, mean)) / apply(nwbl_var, MARGIN = 1, sd, na.rm=T)
SES_wblv

boxplot(SES_wblv, ylim = c(-0.5, 2))
abline(h = 0.0, col = "black", lwd = 3, lty=2)

plot(SES_wblv, pch = 19, cex = 1.5, ylim = c(-0.5, 2))
abline(h = 0.0, col = "black", lwd = 3, lty=2)

hist(as.matrix(nwbl_var), main = "Distribution of expected values - Worker Body Length Variance")
abline(v = cwm.obs.m[,3], col = "blue", lwd = 2)

pval.wblv <- apply(cbind(cwm.obs$wbl_var, nwbl_var), MARGIN = 1, rank)[1,] / 1000
pval.wblv
wblv <- cbind(cwm.obs$wbl_var, nwbl_var)
wblv.m <- as.matrix(t(colMeans(wblv)))
pval.wblv.a <- apply(wblv.m, MARGIN = 1, rank)[1,] / 1000
pval.wblv.a

w.wbl_var <- wilcox.test(SES_wblv, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)#
w.wbl_var


## tongue length - short
SES_tl_0 <- (cwm.obs$tl_0 - apply(ntl_0, MARGIN = 1, mean)) / apply(ntl_0, MARGIN = 1, sd, na.rm=T)
SES_tl_0

boxplot(SES_tl_0)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

plot(SES_tl_0, pch = 19, cex = 1.5)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

hist(as.matrix(ntl_0), main = "Distribution of expected values - Tongue Length")
abline(v = cwm.obs.m[,4], col = "blue", lwd = 2)

pval.tl_0 <- apply(cbind(cwm.obs$tl_0, ntl_0), MARGIN = 1, rank)[1,] / 1000
pval.tl_0
tl_0 <- cbind(cwm.obs$tl_0, ntl_0)
tl_0.m <- as.matrix(t(colMeans(tl_0)))
pval.tl_0.a <- apply(tl_0.m, MARGIN = 1, rank)[1,] / 1000
pval.tl_0.a

w.tl_0 <- wilcox.test(SES_tl_0, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.tl_0


## tongue length - medium
SES_tl_1 <- (cwm.obs$tl_1 - apply(ntl_1, MARGIN = 1, mean)) / apply(ntl_1, MARGIN = 1, sd, na.rm=T)
SES_tl_1

boxplot(SES_tl_1)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

plot(SES_tl_1, pch = 19, cex = 1.5)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

hist(as.matrix(ntl_1), main = "Distribution of expected values - Tongue Length")
abline(v = cwm.obs.m[,5], col = "blue", lwd = 2)

pval.tl_1 <- apply(cbind(cwm.obs$tl_1, ntl_1), MARGIN = 1, rank)[1,] / 1000
pval.tl_1
tl_1 <- cbind(cwm.obs$tl_1, ntl_1)
tl_1.m <- as.matrix(t(colMeans(tl_1)))
pval.tl_1.a <- apply(tl_1.m, MARGIN = 1, rank)[1,] / 1000
pval.tl_1.a

w.tl_1 <- wilcox.test(SES_tl_1, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.tl_1


## tongue length - long
SES_tl_2 <- (cwm.obs$tl_2 - apply(ntl_2, MARGIN = 1, mean)) / apply(ntl_2, MARGIN = 1, sd, na.rm=T)
SES_tl_2

boxplot(SES_tl_2)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

plot(SES_tl_2, pch = 19, cex = 1.5)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

hist(as.matrix(ntl_2), main = "Distribution of expected values - Tongue Length")
abline(v = cwm.obs.m[,6], col = "blue", lwd = 2)

pval.tl_2 <- apply(cbind(cwm.obs$tl_2, ntl_2), MARGIN = 1, rank)[1,] / 1000
pval.tl_2
tl_2 <- cbind(cwm.obs$tl_2, ntl_2)
tl_2.m <- as.matrix(t(colMeans(tl_2)))
pval.tl_2.a <- apply(tl_2.m, MARGIN = 1, rank)[1,] / 1000
pval.tl_2.a

w.tl_2 <- wilcox.test(SES_tl_2, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.tl_2


## nest location - below ground
SES_nest0 <- (cwm.obs$nestl_0 - apply(nnestl_0, MARGIN = 1, mean)) / apply(nnestl_0, MARGIN = 1, sd, na.rm=T)
SES_nest0

boxplot(SES_nest0)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

plot(SES_nest0, pch = 19, cex = 1.5)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

hist(as.matrix(nnestl_0), main = "Distribution of expected values - Belowground Nest Location")
abline(v = cwm.obs.m[,7], col = "blue", lwd = 2)

pval.nest0 <- apply(cbind(cwm.obs$nestl_0, nnestl_0), MARGIN = 1, rank)[1,] / 1000
pval.nest0
nest0 <- cbind(cwm.obs$nestl_0, nnestl_0)
nest0.m <- as.matrix(t(colMeans(nest0)))
pval.nest0.a <- apply(nest0.m, MARGIN = 1, rank)[1,] / 1000
pval.nest0.a

w.nest0 <- wilcox.test(SES_nest0, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.nest0


## nest location - mixed - kleptoparasitic
SES_nest1 <- (cwm.obs$nestl_1 - apply(nnestl_1, MARGIN = 1, mean)) / apply(nnestl_1, MARGIN = 1, sd, na.rm=T)
SES_nest1

boxplot(SES_nest1) #zero not on the y scale!
boxplot(SES_nest1, ylim=c(-1.5, 0.5))
abline(h = 0.0, col = "black", lwd = 3, lty=2)

plot(SES_nest1, pch = 19, cex = 1.5, ylim=c(-1.5, 0.5))
abline(h = 0.0, col = "black", lwd = 3, lty=2)

hist(as.matrix(nnestl_1), main = "Distribution of expected values - Kleptoparasitic Species")
abline(v = cwm.obs.m[,8], col = "blue", lwd = 2)

pval.nest1 <- apply(cbind(cwm.obs$nestl_1, nnestl_1), MARGIN = 1, rank)[1,] / 1000
pval.nest1
nest1 <- cbind(cwm.obs$nestl_1, nnestl_1)
nest1.m <- as.matrix(t(colMeans(nest1)))
pval.nest1.a <- apply(nest1.m, MARGIN = 1, rank)[1,] / 1000
pval.nest1.a

w.nest1 <- wilcox.test(SES_nest1, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)#
w.nest1


## nest location - above ground
SES_nest2 <- (cwm.obs$nestl_2 - apply(nnestl_2, MARGIN = 1, mean)) / apply(nnestl_2, MARGIN = 1, sd, na.rm=T)
SES_nest2

boxplot(SES_nest2)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

plot(SES_nest2, pch = 19, cex = 1.5)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

hist(as.matrix(nnestl_2), main = "Distribution of expected values - Aboveground Nest Location")
abline(v = cwm.obs.m[,9], col = "blue", lwd = 2)

pval.nest2 <- apply(cbind(cwm.obs$nestl_2, nnestl_2), MARGIN = 1, rank)[1,] / 1000
pval.nest2
nest2 <- cbind(cwm.obs$nestl_2, nnestl_2)
nest2.m <- as.matrix(t(colMeans(nest2)))
pval.nest2.a <- apply(nest2.m, MARGIN = 1, rank)[1,] / 1000
pval.nest2.a

w.nest2 <- wilcox.test(SES_nest2, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)#
w.nest2


## intertegular distance
SES_it <- (cwm.obs$it - apply(nit, MARGIN = 1, mean)) / apply(nit, MARGIN = 1, sd, na.rm=T)
SES_it

boxplot(SES_it)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

plot(SES_it, pch = 19, cex = 1.5)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

hist(as.matrix(nit), main = "Distribution of expected values - Intertegular Distance")
abline(v = cwm.obs.m[,10], col = "blue", lwd = 2)

pval.it <- apply(cbind(cwm.obs$it, nit), MARGIN = 1, rank)[1,] / 1000
pval.it
it <- cbind(cwm.obs$it, nit)
it.m <- as.matrix(t(colMeans(it)))
pval.it.a <- apply(it.m, MARGIN = 1, rank)[1,] / 1000
pval.it.a

w.it <- wilcox.test(SES_it, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)#
w.it


## wing length
SES_wingl <- (cwm.obs$rwingl - apply(nwingl, MARGIN = 1, mean)) / apply(nwingl, MARGIN = 1, sd, na.rm=T)
SES_wingl

boxplot(SES_wingl)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

plot(SES_wingl, pch = 19, cex = 1.5)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

hist(as.matrix(nwingl), main = "Distribution of expected values - Wing Length")
abline(v = cwm.obs.m[,11], col = "blue", lwd = 2)

pval.wingl <- apply(cbind(cwm.obs$rwingl, nwingl), MARGIN = 1, rank)[1,] / 1000
pval.wingl
wingl <- cbind(cwm.obs$rwingl, nwingl)
wingl.m <- as.matrix(t(colMeans(wingl)))
pval.wingl.a <- apply(wingl.m, MARGIN = 1, rank)[1,] / 1000
pval.wingl.a

w.wingl <- wilcox.test(SES_wingl, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)#
w.wingl


## head width
SES_headw <- (cwm.obs$rheadw - apply(nheadw, MARGIN = 1, mean)) / apply(nheadw, MARGIN = 1, sd, na.rm=T)
SES_headw

boxplot(SES_headw)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

plot(SES_headw, pch = 19, cex = 1.5)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

hist(as.matrix(nheadw), main = "Distribution of expected values - Head Width")
abline(v = cwm.obs.m[,12], col = "blue", lwd = 2)

pval.headw <- apply(cbind(cwm.obs$rheadw, nheadw), MARGIN = 1, rank)[1,] / 1000
pval.headw
headw <- cbind(cwm.obs$rheadw, nheadw)
headw.m <- as.matrix(t(colMeans(headw)))
pval.headw.a <- apply(headw.m, MARGIN = 1, rank)[1,] / 1000
pval.headw.a

w.headw <- wilcox.test(SES_headw, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)#
w.headw


## eye length
SES_eyel <- (cwm.obs$reyel - apply(neyel, MARGIN = 1, mean)) / apply(neyel, MARGIN = 1, sd, na.rm=T)
SES_eyel

boxplot(SES_eyel)#zero not on the y scale!
boxplot(SES_eyel, ylim=c(-0.5, 2.0))
abline(h = 0.0, col = "black", lwd = 3, lty=2)

plot(SES_eyel, pch = 19, cex = 1.5)
plot(SES_eyel, pch = 19, cex = 1.5, ylim=c(-0.5, 2.0))
abline(h = 0.0, col = "black", lwd = 3, lty=2)

hist(as.matrix(neyel), main = "Distribution of expected values - Eye Length")
abline(v = cwm.obs.m[,13], col = "blue", lwd = 2)

pval.eyel <- apply(cbind(cwm.obs$reyel, neyel), MARGIN = 1, rank)[1,] / 1000
pval.eyel
eyel <- cbind(cwm.obs$reyel, neyel)
eyel.m <- as.matrix(t(colMeans(eyel)))
pval.eyel.a <- apply(eyel.m, MARGIN = 1, rank)[1,] / 1000
pval.eyel.a

w.eyel <- wilcox.test(SES_eyel, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)#
w.eyel


## thorax hair length
SES_thairl <- (cwm.obs$thairl - apply(nthairl, MARGIN = 1, mean)) / apply(nthairl, MARGIN = 1, sd, na.rm=T)
SES_thairl

boxplot(SES_thairl)#zero not on the y scale!
boxplot(SES_thairl, ylim=c(-2.0, 0.5))
abline(h = 0.0, col = "black", lwd = 3, lty=2)

plot(SES_thairl, pch = 19, cex = 1.5)
plot(SES_thairl, pch = 19, cex = 1.5, ylim=c(-2.0, 0.5))
abline(h = 0.0, col = "black", lwd = 3, lty=2)

hist(as.matrix(nthairl), main = "Distribution of expected values - Thorax Hair Length")
abline(v = cwm.obs.m[,14], col = "blue", lwd = 2)

pval.thairl <- apply(cbind(cwm.obs$thairl, nthairl), MARGIN = 1, rank)[1,] / 1000
pval.thairl
thairl <- cbind(cwm.obs$thairl, nthairl)
thairl.m <- as.matrix(t(colMeans(thairl)))
pval.thairl.a <- apply(thairl.m, MARGIN = 1, rank)[1,] / 1000
pval.thairl.a

w.thairl <- wilcox.test(SES_thairl, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)#
w.thairl


## corbicula setae length
SES_setael <- (cwm.obs$bsetael - apply(nbsetael, MARGIN = 1, mean)) / apply(nbsetael, MARGIN = 1, sd, na.rm=T)
SES_setael

boxplot(SES_setael)#zero not on the y scale!
boxplot(SES_setael, ylim=c(-0.5, 1.5))
abline(h = 0.0, col = "black", lwd = 3, lty=2)

plot(SES_setael, pch = 19, cex = 1.5)
plot(SES_setael, pch = 19, cex = 1.5, ylim=c(-0.5, 1.5))
abline(h = 0.0, col = "black", lwd = 3, lty=2)

hist(as.matrix(nbsetael), main = "Distribution of expected values - Corbicula Setae Length")
abline(v = cwm.obs.m[,15], col = "blue", lwd = 2)

pval.setael <- apply(cbind(cwm.obs$bsetael, nbsetael), MARGIN = 1, rank)[1,] / 1000
pval.setael
setael <- cbind(cwm.obs$bsetael, nbsetael)
setael.m <- as.matrix(t(colMeans(setael)))
pval.setael.a <- apply(setael.m, MARGIN = 1, rank)[1,] / 1000
pval.setael.a

w.setael <- wilcox.test(SES_setael, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)#
w.setael


## corbicula length
SES_tibial <- (cwm.obs$rbtl - apply(nbtl, MARGIN = 1, mean)) / apply(nbtl, MARGIN = 1, sd, na.rm=T)
SES_tibial

boxplot(SES_tibial)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

plot(SES_tibial, pch = 19, cex = 1.5)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

hist(as.matrix(nbtl), main = "Distribution of expected values - Corbicula Length")
abline(v = cwm.obs.m[,16], col = "blue", lwd = 2)

pval.tibial <- apply(cbind(cwm.obs$rbtl, nbtl), MARGIN = 1, rank)[1,] / 1000
pval.tibial
tibial <- cbind(cwm.obs$rbtl, nbtl)
tibial.m <- as.matrix(t(colMeans(tibial)))
pval.tibial.a <- apply(tibial.m, MARGIN = 1, rank)[1,] / 1000
pval.tibial.a

w.tibial <- wilcox.test(SES_tibial, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.tibial


## Diversity Indices ##

#Taxonomic diversity

beta.sor <- as.matrix(bb.dist$beta.sor)
beta.sor <- colMeans(beta.sor)

beta.sim <- as.matrix(bb.dist$beta.sim)
beta.sim <- colMeans(beta.sim)

beta.sne <- as.matrix(bb.dist$beta.sne)
beta.sne <- colMeans(beta.sne)

beta.t <- data.frame(beta.sor, beta.sim, beta.sne)
beta.t.m <- as.matrix(t(colMeans(beta.t)))

## taxonomic diversity - beta sor
SES_bsor <- (beta.t$beta.sor - apply(nbsor, MARGIN = 1, mean)) / apply(nbsor, MARGIN = 1, sd, na.rm=T)
SES_bsor

boxplot(SES_bsor)
boxplot(SES_bsor, ylim=c(-14, 0.5))
abline(h = 0.0, col = "black", lwd = 3, lty=2)

plot(SES_bsor, pch = 19, cex = 1.5, ylim=c(-14, 0.5))
abline(h = 0.0, col = "black", lwd = 3, lty=2)

hist(as.matrix(nbsor), main = "Distribution of expected values - Taxonomic Beta-Diversity", xlim=c(0.1, 1))
abline(v = beta.t.m[,1], col = "blue", lwd = 2)

pval.beta.sor <- apply(cbind(beta.t$beta.sor, nbsor), MARGIN = 1, rank)[1,] / 1000
pval.beta.sor
tbsor <- cbind(beta.t$beta.sor, nbsor)
tbsor.m <- as.matrix(t(colMeans(tbsor)))
pval.beta.sor.a <- apply(tbsor.m, MARGIN = 1, rank)[1,] / 1000
pval.beta.sor.a

w.bsor <- wilcox.test(SES_bsor, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.bsor


## taxonomic diversity - beta sim
SES_bsim <- (beta.t$beta.sim - apply(nbsim, MARGIN = 1, mean)) / apply(nbsim, MARGIN = 1, sd, na.rm=T)
SES_bsim

boxplot(SES_bsim)
boxplot(SES_bsim, ylim=c(-14, 0.5))
abline(h = 0.0, col = "black", lwd = 3, lty=2)

plot(SES_bsim, pch = 19, cex = 1.5, ylim=c(-14, 0.5))
abline(h = 0.0, col = "black", lwd = 3, lty=2)

hist(as.matrix(nbsim), main = "Distribution of expected values - Taxonomic Beta-Diversity - Turnover", xlim=c(0.01, 1))
abline(v = beta.t.m[,2], col = "blue", lwd = 2)

pval.beta.sim <- apply(cbind(beta.t$beta.sim, nbsim), MARGIN = 1, rank)[1,] / 1000
pval.beta.sim
tbsim <- cbind(beta.t$beta.sim, nbsim)
tbsim.m <- as.matrix(t(colMeans(tbsim)))
pval.beta.sim.a <- apply(tbsim.m, MARGIN = 1, rank)[1,] / 1000
pval.beta.sim.a

w.bsim <- wilcox.test(SES_bsim, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.bsim

## taxonomic diversity - beta sne
SES_bsne <- (beta.t$beta.sne - apply(nbsne, MARGIN = 1, mean)) / apply(nbsne, MARGIN = 1, sd, na.rm=T)
SES_bsne

boxplot(SES_bsne)
boxplot(SES_bsne, ylim=c(-0.5, 11))
abline(h = 0.0, col = "black", lwd = 3, lty=2)

plot(SES_bsne, pch = 19, cex = 1.5, ylim=c(-0.5, 11))
abline(h = 0.0, col = "black", lwd = 3, lty=2)

hist(as.matrix(nbsne), main = "Distribution of expected values - Taxonomic Beta-Diversity - Nestedness")
abline(v = beta.t.m[,3], col = "blue", lwd = 2)

pval.beta.sne <- apply(cbind(beta.t$beta.sne, nbsne), MARGIN = 1, rank)[1,] / 1000
pval.beta.sne
tbsne <- cbind(beta.t$beta.sne, nbsne)
tbsne.m <- as.matrix(t(colMeans(tbsne)))
pval.beta.sne.a <- apply(tbsne.m, MARGIN = 1, rank)[1,] / 1000
pval.beta.sne.a

w.bsne <- wilcox.test(SES_bsne, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.bsne

# Functional Diversity
fbeta.sor <- as.matrix(bb.fun$funct.beta.sor)
fbeta.sor <- colMeans(fbeta.sor)

fbeta.sim <- as.matrix(bb.fun$funct.beta.sim)
fbeta.sim <- colMeans(fbeta.sim)

fbeta.sne <- as.matrix(bb.fun$funct.beta.sne)
fbeta.sne <- colMeans(fbeta.sne)

beta.f <- data.frame(fbeta.sor, fbeta.sim, fbeta.sne)
beta.f.m <- as.matrix(t(colMeans(beta.f)))

## functional diversity - beta sor
SES_fbsor <- (beta.f$fbeta.sor - apply(nfsor, MARGIN = 1, mean)) / apply(nfsor, MARGIN = 1, sd, na.rm=T)
SES_fbsor

boxplot(SES_fbsor, ylim = c(-5, 0.5))
abline(h = 0.0, col = "black", lwd = 3, lty=2)

plot(SES_fbsor, pch = 19, cex = 1.5, ylim = c(-5, 0.5))
abline(h = 0.0, col = "black", lwd = 3, lty=2)

hist(as.matrix(nfsor), main = "Distribution of expected values - Functional Beta-Diversity", xlim = c(0.1, 1))
abline(v = beta.f.m[,1], col = "blue", lwd = 2)

pval.fbeta.sor <- apply(cbind(beta.f$fbeta.sor, nfsor), MARGIN = 1, rank)[1,] / 1000
pval.fbeta.sor
fbsor <- cbind(beta.f$fbeta.sor, nfsor)
fbsor.m <- as.matrix(t(colMeans(fbsor)))
pval.fbeta.sor.a <- apply(fbsor.m, MARGIN = 1, rank)[1,] / 1000
pval.fbeta.sor.a

w.fbsor <- wilcox.test(SES_fbsor, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.fbsor

## functional diversity - beta sim
SES_fbsim <- (beta.f$fbeta.sim - apply(nfsim, MARGIN = 1, mean)) / apply(nfsim, MARGIN = 1, sd, na.rm=T)
SES_fbsim

boxplot(SES_fbsim)
boxplot(SES_fbsim, ylim = c(-3, 0.5))
abline(h = 0.0, col = "black", lwd = 3, lty=2)

plot(SES_fbsim, pch = 19, cex = 1.5, ylim = c(-3, 0.5))
abline(h = 0.0, col = "black", lwd = 3, lty=2)

hist(as.matrix(nfsim), main = "Distribution of expected values - Functional Beta-Diversity - Turnover", xlim = c(0.01, 1))
abline(v = beta.f.m[,2], col = "blue", lwd = 2)

pval.fbeta.sim <- apply(cbind(beta.f$fbeta.sim, nfsim), MARGIN = 1, rank)[1,] / 1000
pval.fbeta.sim
fbsim <- cbind(beta.f$fbeta.sim, nfsim)
fbsim.m <- as.matrix(t(colMeans(fbsim)))
pval.fbeta.sim.a <- apply(fbsim.m, MARGIN = 1, rank)[1,] / 1000
pval.fbeta.sim.a

w.fbsim <- wilcox.test(SES_fbsim, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.fbsim


## functional diversity - beta sne
SES_fbsne <- (beta.f$fbeta.sne - apply(nfsne, MARGIN = 1, mean)) / apply(nfsne, MARGIN = 1, sd, na.rm=T)
SES_fbsne

boxplot(SES_fbsne)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

plot(SES_fbsne, pch = 19, cex = 1.5)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

hist(as.matrix(nfsne), main = "Distribution of expected values - Functional Beta-Diversity - Nestedness")
abline(v = beta.f.m[,3], col = "blue", lwd = 2)

pval.fbeta.sne <- apply(cbind(beta.f$fbeta.sne, nfsne), MARGIN = 1, rank)[1,] / 1000
pval.fbeta.sne
fbsne <- cbind(beta.f$fbeta.sne, nfsne)
fbsne.m <- as.matrix(t(colMeans(fbsne)))
pval.fbeta.sne.a <- apply(fbsne.m, MARGIN = 1, rank)[1,] / 1000
pval.fbeta.sne.a

w.fbsne <- wilcox.test(SES_fbsne, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.fbsne


# Phylogenetic Diversity
pbeta.sor <- as.matrix(bb.phy$phylo.beta.sor)
pbeta.sor <- colMeans(pbeta.sor)

pbeta.sim <- as.matrix(bb.phy$phylo.beta.sim)
pbeta.sim <- colMeans(pbeta.sim)

pbeta.sne <- as.matrix(bb.phy$phylo.beta.sne)
pbeta.sne <- colMeans(pbeta.sne)

beta.p <- data.frame(pbeta.sor, pbeta.sim, pbeta.sne)
beta.p.m <- as.matrix(t(colMeans(beta.p)))


## phylogenetic diversity - beta sor
SES_pbsor <- (beta.p$pbeta.sor - apply(npsor, MARGIN = 1, mean)) / apply(npsor, MARGIN = 1, sd, na.rm=T)
SES_pbsor

boxplot(SES_pbsor, ylim = c(-8, 0.5))
abline(h = 0.0, col = "black", lwd = 3, lty=2)

plot(SES_pbsor, pch = 19, cex = 1.5, ylim = c(-8, 0.5))
abline(h = 0.0, col = "black", lwd = 3, lty=2)

hist(as.matrix(npsor), main = "Distribution of expected values - Phylogenetic Beta-Diversity", xlim = c(0.1, 1))
abline(v = beta.p.m[,1], col = "blue", lwd = 2)

pval.pbeta.sor <- apply(cbind(beta.p$pbeta.sor, npsor), MARGIN = 1, rank)[1,] / 1000
pval.pbeta.sor
pbsor <- cbind(beta.p$pbeta.sor, npsor)
pbsor.m <- as.matrix(t(colMeans(pbsor)))
pval.pbeta.sor.a <- apply(pbsor.m, MARGIN = 1, rank)[1,] / 1000
pval.pbeta.sor.a

w.pbsor <- wilcox.test(SES_pbsor, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.pbsor

## phylogenetic diversity - beta sim
SES_pbsim <- (beta.p$pbeta.sim - apply(npsim, MARGIN = 1, mean)) / apply(npsim, MARGIN = 1, sd, na.rm=T)
SES_pbsim

boxplot(SES_pbsim)
boxplot(SES_pbsim, ylim = c(-8, 0.5))
abline(h = 0.0, col = "black", lwd = 3, lty=2)

plot(SES_pbsim, pch = 19, cex = 1.5, ylim = c(-8, 0.5))
abline(h = 0.0, col = "black", lwd = 3, lty=2)

hist(as.matrix(npsim), main = "Distribution of expected values - Phylogenetic Beta-Diversity - Turnover", xlim = c(0.01, 1))
abline(v = beta.p.m[,2], col = "blue", lwd = 2)

pval.pbeta.sim <- apply(cbind(beta.p$pbeta.sim, npsim), MARGIN = 1, rank)[1,] / 1000
pval.pbeta.sim
pbsim <- cbind(beta.p$pbeta.sim, npsim)
pbsim.m <- as.matrix(t(colMeans(pbsim)))
pval.pbeta.sim.a <- apply(pbsim.m, MARGIN = 1, rank)[1,] / 1000
pval.pbeta.sim.a

w.pbsim <- wilcox.test(SES_pbsim, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.pbsim

## phylogenetic diversity - beta sne
SES_pbsne <- (beta.p$pbeta.sne - apply(npsne, MARGIN = 1, mean)) / apply(npsne, MARGIN = 1, sd, na.rm=T)
SES_pbsne

boxplot(SES_pbsne)
boxplot(SES_pbsne, ylim = c(-0.5, 5))
abline(h = 0.0, col = "black", lwd = 3, lty=2)

plot(SES_pbsne, pch = 19, cex = 1.5, ylim = c(-0.5, 5))
abline(h = 0.0, col = "black", lwd = 3, lty=2)

hist(as.matrix(npsne), main = "Distribution of expected values - Phylogenetic Beta-Diversity - Nestedness")
abline(v = beta.p.m[,3], col = "blue", lwd = 2)

pval.pbeta.sne <- apply(cbind(beta.p$pbeta.sne, npsne), MARGIN = 1, rank)[1,] / 1000
pval.pbeta.sne
pbsne <- cbind(beta.p$pbeta.sne, npsne)
pbsne.m <- as.matrix(t(colMeans(pbsne)))
pval.pbeta.sne.a <- apply(pbsne.m, MARGIN = 1, rank)[1,] / 1000
pval.pbeta.sne.a

w.pbsne <- wilcox.test(SES_pbsne, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.pbsne

# Phylogenetic Signal

# transpose datasets to calculate SES values
physig <- t(physig)
npsig <- as.data.frame(t(npsig))

# Queen body length variance
SES_physig_qbl_var <-(physig[1,] - apply(npsig[1,], MARGIN = 1, mean)) / apply(npsig[1,], MARGIN = 1, sd, na.rm=T)
SES_physig_qbl_var
pval.physig_qbl_var <- apply(cbind(physig[1,], npsig[1,]), MARGIN = 1, rank)[1,] / 1000
pval.physig_qbl_var
hist(as.matrix(npsig[1,]), main = "Distribution of expected values - Phylogenetic Signal - QBL Variance")
abline(v = physig[1,], col = "blue", lwd = 2)

# Male body length variance
SES_physig_mbl_var <-(physig[2,] - apply(npsig[2,], MARGIN = 1, mean)) / apply(npsig[2,], MARGIN = 1, sd, na.rm=T)
SES_physig_mbl_var
pval.physig_mbl_var <- apply(cbind(physig[2,], npsig[2,]), MARGIN = 1, rank)[1,] / 1000
pval.physig_mbl_var
hist(as.matrix(npsig[2,]), main = "Distribution of expected values - Phylogenetic Signal - MBL Variance")
abline(v = physig[2,], col = "blue", lwd = 2)

# Worker body length variance
SES_physig_wbl_var <-(physig[3,] - apply(npsig[3,], MARGIN = 1, mean)) / apply(npsig[3,], MARGIN = 1, sd, na.rm=T)
SES_physig_wbl_var
pval.physig_wbl_var <- apply(cbind(physig[2,], npsig[2,]), MARGIN = 1, rank)[1,] / 1000
pval.physig_wbl_var
hist(as.matrix(npsig[3,]), main = "Distribution of expected values - Phylogenetic Signal - WBL Variance")
abline(v = physig[3,], col = "blue", lwd = 2)

# Intertegular distance
SES_physig_it <-(physig[4,] - apply(npsig[6,], MARGIN = 1, mean)) / apply(npsig[6,], MARGIN = 1, sd, na.rm=T)
SES_physig_it
pval.physig_it <- apply(cbind(physig[4,], npsig[6,]), MARGIN = 1, rank)[1,] / 1000
pval.physig_it
hist(as.matrix(npsig[6,]), main = "Distribution of expected values - Phylogenetic Signal - IT Distance")
abline(v = physig[4,], col = "blue", lwd = 2)

# Wing length
SES_physig_wl <-(physig[5,] - apply(npsig[7,], MARGIN = 1, mean)) / apply(npsig[7,], MARGIN = 1, sd, na.rm=T)
SES_physig_wl
pval.physig_wl <- apply(cbind(physig[5,], npsig[7,]), MARGIN = 1, rank)[1,] / 1000
pval.physig_wl
hist(as.matrix(npsig[7,]), main = "Distribution of expected values - Phylogenetic Signal - Wing length")
abline(v = physig[5,], col = "blue", lwd = 2)

# Head width
SES_physig_hw <-(physig[6,] - apply(npsig[8,], MARGIN = 1, mean)) / apply(npsig[8,], MARGIN = 1, sd, na.rm=T)
SES_physig_hw
pval.physig_hw <- apply(cbind(physig[6,], npsig[8,]), MARGIN = 1, rank)[1,] / 1000
pval.physig_hw
hist(as.matrix(npsig[8,]), main = "Distribution of expected values - Phylogenetic Signal - Head width")
abline(v = physig[6,], col = "blue", lwd = 2)

# Eye length
SES_physig_el <-(physig[7,] - apply(npsig[9,], MARGIN = 1, mean)) / apply(npsig[9,], MARGIN = 1, sd, na.rm=T)
SES_physig_el
pval.physig_el <- apply(cbind(physig[7,], npsig[9,]), MARGIN = 1, rank)[1,] / 1000
pval.physig_el
hist(as.matrix(npsig[9,]), main = "Distribution of expected values - Phylogenetic Signal - Eye length")
abline(v = physig[7,], col = "blue", lwd = 2)

# Thorax hair length
SES_physig_hl <-(physig[8,] - apply(npsig[10,], MARGIN = 1, mean)) / apply(npsig[10,], MARGIN = 1, sd, na.rm=T)
SES_physig_hl
pval.physig_hl <- apply(cbind(physig[8,], npsig[10,]), MARGIN = 1, rank)[1,] / 1000
pval.physig_hl
hist(as.matrix(npsig[10,]), main = "Distribution of expected values - Phylogenetic Signal - Hair length")
abline(v = physig[8,], col = "blue", lwd = 2)

# Corbicula setae length
SES_physig_sl <-(physig[9,] - apply(npsig[11,], MARGIN = 1, mean)) / apply(npsig[11,], MARGIN = 1, sd, na.rm=T)
SES_physig_sl
pval.physig_sl <- apply(cbind(physig[9,], npsig[11,]), MARGIN = 1, rank)[1,] / 1000
pval.physig_sl
hist(as.matrix(npsig[11,]), main = "Distribution of expected values - Phylogenetic Signal - Setae length")
abline(v = physig[9,], col = "blue", lwd = 2)

# Corbicula length
SES_physig_cl <-(physig[10,] - apply(npsig[12,], MARGIN = 1, mean)) / apply(npsig[12,], MARGIN = 1, sd, na.rm=T)
SES_physig_cl
pval.physig_cl <- apply(cbind(physig[10,], npsig[12,]), MARGIN = 1, rank)[1,] / 1000
pval.physig_cl
hist(as.matrix(npsig[12,]), main = "Distribution of expected values - Phylogenetic Signal - Corbicula length")
abline(v = physig[10,], col = "blue", lwd = 2)


###################################################################################
####Differences in Overall Observed SES from Null Communities######################

SES <- cbind(SES_qblv, SES_mblv, SES_wblv, SES_tl_0, SES_tl_1, SES_tl_2, SES_nest0, SES_nest1,
             SES_nest2, SES_it, SES_wingl, SES_headw, SES_eyel, SES_thairl, SES_setael, SES_tibial,
             SES_bsor, SES_bsim, SES_bsne, SES_fbsor, SES_fbsim, SES_fbsne, SES_pbsor, SES_pbsim, SES_pbsne)
SES <- as.data.frame(SES)
str(SES)


write.csv(SES, file = "SES_Regional.csv")

# load the dataset
SES <- read.csv("SES_Regional.csv", row.names=1)

###################################################################################
##Differences in Observed SES from Null Communities by Treatment###################

### Group by trait categories

library(viridis)
library(reshape2)

## Body Size Traits ###############################################################
SES_bl <- as.data.frame(SES[,c(3,2,1,10)])
colnames(SES_bl) <- c("Worker BL Variance", "Male BL Variance", "Queen BL Variance", "Inter-tegular Distance")
SES_bl$trmt <- c("U", "A", "U", "U", "U", "A", "A", "U", "A", "U", "U", "U", "U", "A", "U",
              "A", "U", "U", "A", "A")
SES_bl$trmt <- as.factor(SES_bl$trmt)
str(SES_bl)

urban.bl <- SES_bl[which(SES_bl$trmt == "U"),]
ag.bl <- SES_bl[which(SES_bl$trmt == "A"),]

urban.bl.SES <- melt(urban.bl)
colnames(urban.bl.SES) <- c("trmt","trait","ses")

ag.bl.SES <- melt(ag.bl)
colnames(ag.bl.SES) <- c("trmt","trait","ses")


#urban
par(mar = c(5,9,4,2))
boxplot(ses ~ trait, data = urban.bl.SES, col = viridis(4, alpha = 0.6),
        xlab = "Standardized Effect Sizes (SES)", ylab = "",
        horizontal = TRUE, las = 1, range = 0, main = "Urban")
stripchart(ses ~ trait, data = urban.bl.SES, col = viridis(4),
           pch = 19, cex = 2, las = 1, add = TRUE, method = "jitter", jitter = 0.2)
abline(v = 0.0, col = "black", lwd = 3, lty=2)

w.wbl_var.u <- wilcox.test(urban.bl[,1], y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)#
w.wbl_var.u

w.mbl_var.u <- wilcox.test(urban.bl[,2], y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)#
w.mbl_var.u

w.qbl_var.u <- wilcox.test(urban.bl[,3], y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)#
w.qbl_var.u

w.it.u <- wilcox.test(urban.bl[,4], y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)#
w.it.u

#ag
par(mar = c(5,9,4,2))
boxplot(ses ~ trait, data = ag.bl.SES, col = viridis(4, alpha = 0.6),
        xlab = "Standardized Effect Sizes (SES)", ylab = "",
        horizontal = TRUE, las = 1, range = 0, main = "Agriculture")
stripchart(ses ~ trait, data = ag.bl.SES, col = viridis(4),
           pch = 19, cex = 2, las = 1, add = TRUE, method = "jitter", jitter = 0.2)
abline(v = 0.0, col = "black", lwd = 3, lty=2)

w.wbl_var.a <- wilcox.test(ag.bl[,1], y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)#
w.wbl_var.a

w.mbl_var.a <- wilcox.test(ag.bl[,2], y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)#
w.mbl_var.a

w.qbl_var.a <- wilcox.test(ag.bl[,3], y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)#
w.qbl_var.a

w.it.a <- wilcox.test(ag.bl[,4], y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)#
w.it.a

## Tongue Length ###############################################################
SES_tl <- as.data.frame(SES[,c(6,5,4)])
colnames(SES_tl) <- c("Long Tongue", "Medium Tongue", "Short Tongue")
SES_tl$trmt <- c("U", "A", "U", "U", "U", "A", "A", "U", "A", "U", "U", "U", "U", "A", "U",
                 "A", "U", "U", "A", "A")
SES_tl$trmt <- as.factor(SES_tl$trmt)
str(SES_tl)

urban.tl <- SES_tl[which(SES_tl$trmt == "U"),]
ag.tl <- SES_tl[which(SES_tl$trmt == "A"),]

urban.tl.SES <- melt(urban.tl)
colnames(urban.tl.SES) <- c("trmt","trait","ses")

ag.tl.SES <- melt(ag.tl)
colnames(ag.tl.SES) <- c("trmt","trait","ses")


#urban
par(mar = c(5,9,4,2))
boxplot(ses ~ trait, data = urban.tl.SES, col = viridis(3, alpha = 0.6),
        xlab = "Standardized Effect Sizes (SES)", ylab = "",
        horizontal = TRUE, las = 1, range = 0, main = "Urban")
stripchart(ses ~ trait, data = urban.tl.SES, col = viridis(3),
           pch = 19, cex = 2, las = 1, add = TRUE, method = "jitter", jitter = 0.2)
abline(v = 0.0, col = "black", lwd = 3, lty=2)

w.tl0.u <- wilcox.test(urban.tl[,3], y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.tl0.u

w.tl1.u <- wilcox.test(urban.tl[,2], y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.tl1.u

w.tl2.u <- wilcox.test(urban.tl[,1], y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.tl2.u

#ag
par(mar = c(5,9,4,2))
boxplot(ses ~ trait, data = ag.tl.SES, col = viridis(3, alpha = 0.6),
        xlab = "Standardized Effect Sizes (SES)", ylab = "",
        horizontal = TRUE, las = 1, range = 0, main = "Agriculture")
stripchart(ses ~ trait, data = ag.tl.SES, col = viridis(3),
           pch = 19, cex = 2, las = 1, add = TRUE, method = "jitter", jitter = 0.2)
abline(v = 0.0, col = "black", lwd = 3, lty=2)

w.tl0.a <- wilcox.test(ag.tl[,3], y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.tl0.a

w.tl1.a <- wilcox.test(ag.tl[,2], y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.tl1.a

w.tl2.a <- wilcox.test(ag.tl[,1], y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.tl2.a


## Nest Location ###############################################################
SES_nest <- as.data.frame(SES[,c(8,7,9)])
colnames(SES_nest) <- c("Kleptoparasitic", "Belowground Nests", "Aboveground Nests")
SES_nest$trmt <- c("U", "A", "U", "U", "U", "A", "A", "U", "A", "U", "U", "U", "U", "A", "U",
                 "A", "U", "U", "A", "A")
SES_nest$trmt <- as.factor(SES_nest$trmt)
str(SES_nest)

urban.nest <- SES_nest[which(SES_nest$trmt == "U"),]
ag.nest <- SES_nest[which(SES_nest$trmt == "A"),]

urban.nest.SES <- melt(urban.nest)
colnames(urban.nest.SES) <- c("trmt","trait","ses")

ag.nest.SES <- melt(ag.nest)
colnames(ag.nest.SES) <- c("trmt","trait","ses")


#urban
par(mar = c(5,9,4,2))
boxplot(ses ~ trait, data = urban.nest.SES, col = viridis(3, alpha = 0.6),
        xlab = "Standardized Effect Sizes (SES)", ylab = "",
        horizontal = TRUE, las = 1, range = 0, main = "Urban")
stripchart(ses ~ trait, data = urban.nest.SES, col = viridis(3),
           pch = 19, cex = 2, las = 1, add = TRUE, method = "jitter", jitter = 0.2)
abline(v = 0.0, col = "black", lwd = 3, lty=2)

w.nest0.u <- wilcox.test(urban.nest[,2], y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.nest0.u

w.nest1.u <- wilcox.test(urban.nest[,1], y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.nest1.u

w.nest2.u <- wilcox.test(urban.nest[,3], y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.nest2.u

#ag
par(mar = c(5,9,4,2))
boxplot(ses ~ trait, data = ag.nest.SES, col = viridis(3, alpha = 0.6),
        xlab = "Standardized Effect Sizes (SES)", ylab = "",
        horizontal = TRUE, las = 1, range = 0, main = "Agriculture")
stripchart(ses ~ trait, data = ag.nest.SES, col = viridis(3),
           pch = 19, cex = 2, las = 1, add = TRUE, method = "jitter", jitter = 0.2)
abline(v = 0.0, col = "black", lwd = 3, lty=2)

w.nest0.a <- wilcox.test(ag.nest[,2], y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.nest0.a

w.nest1.a <- wilcox.test(ag.nest[,1], y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.nest1.a

w.nest2.a <- wilcox.test(ag.nest[,3], y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.nest2.a


## Other Traits ###############################################################
SES_trait <- as.data.frame(SES[,c(15,16,14,13,12,11)])
colnames(SES_trait) <- c("Setae Length", "Tibia Length", "Thorax Hair Length", "Eye Length", "Head Width", "Wing Length")
SES_trait$trmt <- c("U", "A", "U", "U", "U", "A", "A", "U", "A", "U", "U", "U", "U", "A", "U",
                   "A", "U", "U", "A", "A")
SES_trait$trmt <- as.factor(SES_trait$trmt)
str(SES_trait)

urban.trait <- SES_trait[which(SES_trait$trmt == "U"),]
ag.trait <- SES_trait[which(SES_trait$trmt == "A"),]

urban.trait.SES <- melt(urban.trait)
colnames(urban.trait.SES) <- c("trmt","trait","ses")

ag.trait.SES <- melt(ag.trait)
colnames(ag.trait.SES) <- c("trmt","trait","ses")


#urban
par(mar = c(5,9,4,2))
boxplot(ses ~ trait, data = urban.trait.SES, col = viridis(6, alpha = 0.6),
        xlab = "Standardized Effect Sizes (SES)", ylab = "", ylim = c(-2,2),
        horizontal = TRUE, las = 1, range = 0, main = "Urban")
stripchart(ses ~ trait, data = urban.trait.SES, col = viridis(6),
           pch = 19, cex = 2, las = 1, add = TRUE, method = "jitter", jitter = 0.2)
abline(v = 0.0, col = "black", lwd = 3, lty=2)

w.wingl.u <- wilcox.test(urban.trait[,6], y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)#
w.wingl.u

w.headw.u <- wilcox.test(urban.trait[,5], y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)#
w.headw.u

w.eyel.u <- wilcox.test(urban.trait[,4], y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)#
w.eyel.u

w.thairl.u <- wilcox.test(urban.trait[,3], y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)#
w.thairl.u

w.setael.u <- wilcox.test(urban.trait[,1], y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)#
w.setael.u

w.tibial.u <- wilcox.test(urban.trait[,2], y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.tibial.u

#ag
par(mar = c(5,9,4,2))
boxplot(ses ~ trait, data = ag.trait.SES, col = viridis(6, alpha = 0.6),
        xlab = "Standardized Effect Sizes (SES)", ylab = "", ylim = c(-2,2),
        horizontal = TRUE, las = 1, range = 0, main = "Agriculture")
stripchart(ses ~ trait, data = ag.trait.SES, col = viridis(6),
           pch = 19, cex = 2, las = 1, add = TRUE, method = "jitter", jitter = 0.2)
abline(v = 0.0, col = "black", lwd = 3, lty=2)

w.wingl.a <- wilcox.test(ag.trait[,6], y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)#
w.wingl.a

w.headw.a <- wilcox.test(ag.trait[,5], y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)#
w.headw.a

w.eyel.a <- wilcox.test(ag.trait[,4], y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)#
w.eyel.a

w.thairl.a <- wilcox.test(ag.trait[,3], y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)#
w.thairl.a

w.setael.a <- wilcox.test(ag.trait[,1], y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)#
w.setael.a

w.tibial.a <- wilcox.test(ag.trait[,2], y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.tibial.a


## Taxonomic Beta-diversity ###############################################################
SES_tbeta <- as.data.frame(SES[,c(19,18,17)])
colnames(SES_tbeta) <- c("Nestedness", "Turnover", "Total")
SES_tbeta$trmt <- c("U", "A", "U", "U", "U", "A", "A", "U", "A", "U", "U", "U", "U", "A", "U",
                   "A", "U", "U", "A", "A")
SES_tbeta$trmt <- as.factor(SES_tbeta$trmt)
str(SES_tbeta)

urban.tbeta <- SES_tbeta[which(SES_tbeta$trmt == "U"),]
ag.tbeta <- SES_tbeta[which(SES_tbeta$trmt == "A"),]

urban.tbeta.SES <- melt(urban.tbeta)
colnames(urban.tbeta.SES) <- c("trmt","trait","ses")

ag.tbeta.SES <- melt(ag.tbeta)
colnames(ag.tbeta.SES) <- c("trmt","trait","ses")


#urban
par(mar = c(5,9,4,2))
boxplot(ses ~ trait, data = urban.tbeta.SES, col = viridis(3, alpha = 0.6),
        xlab = "Standardized Effect Sizes (SES)", ylab = "", ylim = c(-15,15),
        horizontal = TRUE, las = 1, range = 0, main = "Urban")
stripchart(ses ~ trait, data = urban.tbeta.SES, col = viridis(3),
           pch = 19, cex = 2, las = 1, add = TRUE, method = "jitter", jitter = 0.2)
abline(v = 0.0, col = "black", lwd = 3, lty=2)

w.bsor.u <- wilcox.test(urban.tbeta[,3], y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.bsor.u

w.bsim.u <- wilcox.test(urban.tbeta[,2], y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.bsim.u

w.bsne.u <- wilcox.test(urban.tbeta[,1], y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.bsne.u

#ag
par(mar = c(5,9,4,2))
boxplot(ses ~ trait, data = ag.tbeta.SES, col = viridis(3, alpha = 0.6),
        xlab = "Standardized Effect Sizes (SES)", ylab = "", ylim = c(-15,15),
        horizontal = TRUE, las = 1, range = 0, main = "Agriculture")
stripchart(ses ~ trait, data = ag.tbeta.SES, col = viridis(3),
           pch = 19, cex = 2, las = 1, add = TRUE, method = "jitter", jitter = 0.2)
abline(v = 0.0, col = "black", lwd = 3, lty=2)

w.bsor.a <- wilcox.test(ag.tbeta[,3], y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.bsor.a

w.bsim.a <- wilcox.test(ag.tbeta[,2], y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.bsim.a

w.bsne.a <- wilcox.test(ag.tbeta[,1], y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.bsne.a


## Functional Beta-diversity ###############################################################
SES_fbeta <- as.data.frame(SES[,c(22,21,20)])
colnames(SES_fbeta) <- c("Nestedness", "Turnover", "Total")
SES_fbeta$trmt <- c("U", "A", "U", "U", "U", "A", "A", "U", "A", "U", "U", "U", "U", "A", "U",
                    "A", "U", "U", "A", "A")
SES_fbeta$trmt <- as.factor(SES_fbeta$trmt)
str(SES_fbeta)

urban.fbeta <- SES_fbeta[which(SES_fbeta$trmt == "U"),]
ag.fbeta <- SES_fbeta[which(SES_fbeta$trmt == "A"),]

urban.fbeta.SES <- melt(urban.fbeta)
colnames(urban.fbeta.SES) <- c("trmt","trait","ses")

ag.fbeta.SES <- melt(ag.fbeta)
colnames(ag.fbeta.SES) <- c("trmt","trait","ses")


#urban
par(mar = c(5,9,4,2))
boxplot(ses ~ trait, data = urban.fbeta.SES, col = viridis(3, alpha = 0.6),
        xlab = "Standardized Effect Sizes (SES)", ylab = "", ylim = c(-6, 6),
        horizontal = TRUE, las = 1, range = 0, main = "Urban")
stripchart(ses ~ trait, data = urban.fbeta.SES, col = viridis(3),
           pch = 19, cex = 2, las = 1, add = TRUE, method = "jitter", jitter = 0.2)
abline(v = 0.0, col = "black", lwd = 3, lty=2)

w.fbsor.u <- wilcox.test(urban.fbeta[,3], y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.fbsor.u

w.fbsim.u <- wilcox.test(urban.fbeta[,2], y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.fbsim.u

w.fbsne.u <- wilcox.test(urban.fbeta[,1], y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.fbsne.u

#ag
par(mar = c(5,9,4,2))
boxplot(ses ~ trait, data = ag.fbeta.SES, col = viridis(3, alpha = 0.6),
        xlab = "Standardized Effect Sizes (SES)", ylab = "", ylim = c(-6, 6),
        horizontal = TRUE, las = 1, range = 0, main = "Agriculture")
stripchart(ses ~ trait, data = ag.fbeta.SES, col = viridis(3),
           pch = 19, cex = 2, las = 1, add = TRUE, method = "jitter", jitter = 0.2)
abline(v = 0.0, col = "black", lwd = 3, lty=2)

w.fbsor.a <- wilcox.test(ag.fbeta[,3], y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.fbsor.a

w.fbsim.a <- wilcox.test(ag.fbeta[,2], y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.fbsim.a

w.fbsne.a <- wilcox.test(ag.fbeta[,1], y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.fbsne.a


## Phylogenetic Beta-diversity ###############################################################
SES_pbeta <- as.data.frame(SES[,c(25,24,23)])
colnames(SES_pbeta) <- c("Nestedness", "Turnover", "Total")
SES_pbeta$trmt <- c("U", "A", "U", "U", "U", "A", "A", "U", "A", "U", "U", "U", "U", "A", "U",
                    "A", "U", "U", "A", "A")
SES_pbeta$trmt <- as.factor(SES_pbeta$trmt)
str(SES_pbeta)

urban.pbeta <- SES_pbeta[which(SES_pbeta$trmt == "U"),]
ag.pbeta <- SES_pbeta[which(SES_pbeta$trmt == "A"),]

urban.pbeta.SES <- melt(urban.pbeta)
colnames(urban.pbeta.SES) <- c("trmt","trait","ses")

ag.pbeta.SES <- melt(ag.pbeta)
colnames(ag.pbeta.SES) <- c("trmt","trait","ses")

#urban
par(mar = c(5,9,4,2))
boxplot(ses ~ trait, data = urban.pbeta.SES, col = viridis(3, alpha = 0.6),
        xlab = "Standardized Effect Sizes (SES)", ylab = "", ylim = c(-8, 5),
        horizontal = TRUE, las = 1, range = 0, main = "Urban")
stripchart(ses ~ trait, data = urban.pbeta.SES, col = viridis(3),
           pch = 19, cex = 2, las = 1, add = TRUE, method = "jitter", jitter = 0.2)
abline(v = 0.0, col = "black", lwd = 3, lty=2)

w.pbsor.u <- wilcox.test(urban.pbeta[,3], y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.pbsor.u

w.pbsim.u <- wilcox.test(urban.pbeta[,2], y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.pbsim.u

w.pbsne.u <- wilcox.test(urban.pbeta[,1], y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.pbsne.u

#ag
par(mar = c(5,9,4,2))
boxplot(ses ~ trait, data = ag.pbeta.SES, col = viridis(3, alpha = 0.6),
        xlab = "Standardized Effect Sizes (SES)", ylab = "", ylim = c(-8, 5),
        horizontal = TRUE, las = 1, range = 0, main = "Agriculture")
stripchart(ses ~ trait, data = ag.pbeta.SES, col = viridis(3),
           pch = 19, cex = 2, las = 1, add = TRUE, method = "jitter", jitter = 0.2)
abline(v = 0.0, col = "black", lwd = 3, lty=2)

w.pbsor.a <- wilcox.test(ag.pbeta[,3], y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.pbsor.a

w.pbsim.a <- wilcox.test(ag.pbeta[,2], y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.pbsim.a

w.pbsne.a <- wilcox.test(ag.pbeta[,1], y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.pbsne.a


#########################################################################################################

### figure panel - all diversity indices

library(ggplot2)
library(ggthemes)

SES_tbeta$metric <- rep(c("Taxonomic"),each = 20)
SES_fbeta$metric <- rep(c("Functional"),each = 20)
SES_pbeta$metric <- rep(c("Phylogenetic"),each = 20)

SES_div <- rbind(SES_tbeta, SES_fbeta, SES_pbeta)
SES_div.m <- melt(SES_div)
colnames(SES_div.m) <- c("trmt","metric", "var", "ses")

levels(SES_div.m$trmt)[levels(SES_div.m$trmt)=="U"] <- "Urban"
levels(SES_div.m$trmt)[levels(SES_div.m$trmt)=="A"] <- "Agriculture"

SES_div.m$metric <- factor(SES_div.m$metric, levels = c("Taxonomic", "Functional", "Phylogenetic"))


png("SES_Div.png", width = 2000, height = 1000, pointsize = 20)

ggplot(SES_div.m, aes(x=ses, y=var, fill = var)) +
              geom_boxplot(outlier.shape = NA) +
              geom_dotplot(position = position_jitter(width = 0.2, height = 0.2),
                           dotsize = 2,
                           binaxis = "y",
                           stackdir = "center") +
              facet_grid(metric ~ trmt) + 
              coord_cartesian(xlim = c(-12, 12)) +
              theme_few() +
              theme(text = element_text(size = 24, color = "black"),
                    axis.text.x = element_text(color = "black"),
                    axis.text.y = element_text(color = "black"),
                    axis.title.x = element_text(vjust = 1, size = 28),
                    axis.title.y = element_text(vjust = 1, size = 28),
                    strip.text = element_text(face = "bold", size = rel(1.1)),
                    strip.background = element_rect(fill = "gray88", color = "black"),
                    legend.position = "none") +
              labs(x = "Standardized Effect Sizes (SES)", y = "Beta-diversity Metrics") +
              geom_vline(xintercept = 0, size = 1.2, linetype = "dashed") +
              scale_fill_viridis(alpha = 0.7, discrete = TRUE, option = "D")

dev.off()

#########################################################################################
## Figure panel version 2

# load the dataset
SES <- read.csv("SES_Regional.csv", row.names=1)

SES$trmt <- c("Urban", "Agriculture", "Urban", "Urban", "Urban", "Agriculture", "Agriculture",
              "Urban", "Agriculture", "Urban", "Urban", "Urban", "Urban", "Agriculture", "Urban",
              "Agriculture", "Urban", "Urban", "Agriculture", "Agriculture")
str(SES)
SES$trmt <- as.factor(SES$trmt)
str(SES)

# All three together
png("Regional.SES.png", width = 2800, height = 3000, pointsize = 30)

par(mfrow=c(3,3))
par(mar=c(5,8,4,2))

#Taxonomic
boxplot(SES$SES_bsor ~ trmt, data = SES, col = c("gray60", "orchid"),
        ylab = "Standardized Effect Sizes (SES)", xlab = "", cex.main = 2,
        cex.lab = 2.2, cex.axis = 2, ylim = c(-12,12), outline = FALSE)
stripchart(SES$SES_bsor ~ trmt, data = SES, col = c("gray28", "orchid4"), vertical = TRUE,
           pch = 19, cex = 2, add = TRUE, method = "jitter", jitter = 0.2)
abline(h = 0.0, col = "black", lwd = 3, lty=2)
text(0.8, 11.5, "Total", pos = 2, font = 2, cex = 2)
text(1.1, -7, "**", pos = 2, font = 2, cex = 3)
text(2.13, -7, "***", pos = 2, font = 2, cex = 3)

boxplot(SES$SES_bsim ~ trmt, data = SES, col = c("gray60", "orchid"),
        ylab = "Standardized Effect Sizes (SES)", xlab = "", main = "Taxonomic", cex.main = 2.5,
        cex.lab = 2.2, cex.axis = 2, ylim = c(-12,12), outline = FALSE)
stripchart(SES$SES_bsim ~ trmt, data = SES, col = c("gray28", "orchid4"), vertical = TRUE,
           pch = 19, cex = 2, add = TRUE, method = "jitter", jitter = 0.2)
abline(h = 0.0, col = "black", lwd = 3, lty=2)
text(0.73, 10.8, "Turnover", pos = 3, font = 2, cex = 2)
text(1.1, -7, "**", pos = 2, font = 2, cex = 3)
text(2.13, -7, "***", pos = 2, font = 2, cex = 3)

boxplot(SES$SES_bsne ~ trmt, data = SES, col = c("gray60", "orchid"),
        ylab = "Standardized Effect Sizes (SES)", xlab = "", cex.main = 2,
        cex.lab = 2.2, cex.axis = 2, ylim = c(-12,12), outline = FALSE)
stripchart(SES$SES_bsne ~ trmt, data = SES, col = c("gray28", "orchid4"), vertical = TRUE,
           pch = 19, cex = 2, add = TRUE, method = "jitter", jitter = 0.2)
abline(h = 0.0, col = "black", lwd = 3, lty=2)
text(0.82, 10.8, "Nestedness", pos = 3, font = 2, cex = 2)
text(1.1, 9.5, "**", pos = 2, font = 2, cex = 3)
text(2.13, 9.5, "***", pos = 2, font = 2, cex = 3)

#Phylogenetic
boxplot(SES$SES_pbsor ~ trmt, data = SES, col = c("gray60", "orchid"),
        ylab = "Standardized Effect Sizes (SES)", xlab = "", cex.main = 2,
        cex.lab = 2.2, cex.axis = 2, ylim = c(-12,12), outline = FALSE)
stripchart(SES$SES_pbsor ~ trmt, data = SES, col = c("gray28", "orchid4"), vertical = TRUE,
           pch = 19, cex = 2, add = TRUE, method = "jitter", jitter = 0.2)
abline(h = 0.0, col = "black", lwd = 3, lty=2)
text(0.8, 11.5, "Total", pos = 2, font = 2, cex = 2)
text(1.1, -2, "**", pos = 2, font = 2, cex = 3)
text(2.13, -2, "***", pos = 2, font = 2, cex = 3)

boxplot(SES$SES_pbsim ~ trmt, data = SES, col = c("gray60", "orchid"),
        ylab = "Standardized Effect Sizes (SES)", xlab = "", main = "Phylogenetic", cex.main = 2.5,
        cex.lab = 2.2, cex.axis = 2, ylim = c(-12,12), outline = FALSE)
stripchart(SES$SES_pbsim ~ trmt, data = SES, col = c("gray28", "orchid4"), vertical = TRUE,
           pch = 19, cex = 2, add = TRUE, method = "jitter", jitter = 0.2)
abline(h = 0.0, col = "black", lwd = 3, lty=2)
text(0.73, 10.8, "Turnover", pos = 3, font = 2, cex = 2)
text(1.1, -2, "**", pos = 2, font = 2, cex = 3)
text(2.13, -2, "***", pos = 2, font = 2, cex = 3)

boxplot(SES$SES_pbsne ~ trmt, data = SES, col = c("gray60", "orchid"),
        ylab = "Standardized Effect Sizes (SES)", xlab = "", cex.main = 2,
        cex.lab = 2.2, cex.axis = 2, ylim = c(-12,12), outline = FALSE)
stripchart(SES$SES_pbsne ~ trmt, data = SES, col = c("gray28", "orchid4"), vertical = TRUE,
           pch = 19, cex = 2, add = TRUE, method = "jitter", jitter = 0.2)
abline(h = 0.0, col = "black", lwd = 3, lty=2)
text(0.82, 10.8, "Nestedness", pos = 3, font = 2, cex = 2)
text(1.1, 5.5, "**", pos = 2, font = 2, cex = 3)
text(2.13, 5.5, "***", pos = 2, font = 2, cex = 3)

#Functional
boxplot(SES$SES_fbsor ~ trmt, data = SES, col = c("gray60", "orchid"),
        ylab = "Standardized Effect Sizes (SES)", xlab = "", cex.main = 2,
        cex.lab = 2.2, cex.axis = 2, ylim = c(-12,12), outline = FALSE)
stripchart(SES$SES_fbsor ~ trmt, data = SES, col = c("gray28", "orchid4"), vertical = TRUE,
           pch = 19, cex = 2, add = TRUE, method = "jitter", jitter = 0.2)
abline(h = 0.0, col = "black", lwd = 3, lty=2)
text(0.8, 11.5, "Total", pos = 2, font = 2, cex = 2)
text(1.1, 1, "**", pos = 2, font = 2, cex = 3)
text(2.13, -1.5, "***", pos = 2, font = 2, cex = 3)

boxplot(SES$SES_fbsim ~ trmt, data = SES, col = c("gray60", "orchid"),
        ylab = "Standardized Effect Sizes (SES)", xlab = "", main = "Functional", cex.main = 2.5,
        cex.lab = 2.2, cex.axis = 2, ylim = c(-12,12), outline = FALSE)
stripchart(SES$SES_fbsim ~ trmt, data = SES, col = c("gray28", "orchid4"), vertical = TRUE,
           pch = 19, cex = 2, add = TRUE, method = "jitter", jitter = 0.2)
abline(h = 0.0, col = "black", lwd = 3, lty=2)
text(0.73, 10.8, "Turnover", pos = 3, font = 2, cex = 2)
text(1.1, 0.5, "**", pos = 2, font = 2, cex = 3)
text(2.13, 0.5, "***", pos = 2, font = 2, cex = 3)

boxplot(SES$SES_fbsne ~ trmt, data = SES, col = c("gray60", "orchid"),
        ylab = "Standardized Effect Sizes (SES)", xlab = "", cex.main = 2,
        cex.lab = 2.2, cex.axis = 2, ylim = c(-12,12), outline = FALSE)
stripchart(SES$SES_fbsne ~ trmt, data = SES, col = c("gray28", "orchid4"), vertical = TRUE,
           pch = 19, cex = 2, add = TRUE, method = "jitter", jitter = 0.2)
abline(h = 0.0, col = "black", lwd = 3, lty=2)
text(0.82, 10.8, "Nestedness", pos = 3, font = 2, cex = 2)
text(2.13, 1, "***", pos = 2, font = 2, cex = 3)

dev.off()


### figure panel - all CWMs

SES_cwm <- cbind(SES_qblv, SES_mblv, SES_wblv, SES_it, SES_tl_0, SES_tl_1, SES_tl_2, SES_nest0, SES_nest1,
             SES_nest2, SES_wingl, SES_headw, SES_eyel, SES_thairl, SES_setael, SES_tibial)
SES_cwm <- as.data.frame(SES_cwm)
str(SES_cwm)
SES_cwm$trmt <- c("U", "A", "U", "U", "U", "A", "A", "U", "A", "U", "U", "U", "U", "A", "U",
                    "A", "U", "U", "A", "A")
SES_cwm$trmt <- as.factor(SES_cwm$trmt)
str(SES_cwm)

SES_cwm.m <- melt(SES_cwm)
colnames(SES_cwm.m) <- c("trmt","cwm", "ses")

levels(SES_cwm.m$trmt)[levels(SES_cwm.m$trmt)=="U"] <- "Urban"
levels(SES_cwm.m$trmt)[levels(SES_cwm.m$trmt)=="A"] <- "Agriculture"

str(SES_cwm.m)

png("SES_CWM.png", width = 2000, height = 1200, pointsize = 20)

ggplot(SES_cwm.m, aes(x=ses, y=cwm, fill = cwm)) +
  geom_boxplot(outlier.shape = NA) +
  geom_dotplot(position = position_jitter(width = 0.2, height = 0.2),
               dotsize = 0.6,
               binaxis = "y",
               stackdir = "center") +
  facet_grid( ~ trmt) + 
  coord_cartesian(xlim = c(-3, 3)) +
  theme_few() +
  theme(text = element_text(size = 24, color = "black"),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.title.x = element_text(vjust = 1, size = 28),
        axis.title.y = element_text(vjust = 1, size = 28),
        strip.text = element_text(face = "bold", size = rel(1.1)),
        strip.background = element_rect(fill = "gray88", color = "black"),
        legend.position = "none") +
  labs(x = "Standardized Effect Sizes (SES)", y = "Community-weighted Means") +
  geom_vline(xintercept = 0, size = 1.2, linetype = "dashed") +
  scale_fill_viridis(alpha = 0.7, discrete = TRUE, option = "D")

  dev.off()
  
#########################################################################################################
  
### figure panel - Total diversity indices
  
SES$trmt <- c("U", "A", "U", "U", "U", "A", "A", "U", "A", "U", "U", "U", "U", "A", "U",
                    "A", "U", "U", "A", "A")
  
levels(SES$trmt)[levels(SES$trmt)=="U"] <- "Urban"
levels(SES$trmt)[levels(SES$trmt)=="A"] <- "Agriculture"

SES$trmt <- as.factor(SES$trmt)
str(SES)

# Total beta-diversity
png("R_SES_TotalDiv.png", width = 2500, height = 1000, pointsize = 30)

par(mfrow=c(1,3))
par(mar=c(5,8,4,2))

boxplot(SES$SES_bsor ~ trmt, data = SES, col = c("gray60", "aquamarine3"),
        ylab = "Standardized Effect Sizes (SES)", xlab = "", main = "Taxonomic", cex.main = 2,
        cex.lab = 1.6, cex.axis = 1.2, ylim = c(-12,5), outline = FALSE)
stripchart(SES$SES_bsor ~ trmt, data = SES, col = c("gray28", "aquamarine4"), vertical = TRUE,
           pch = 19, cex = 2, add = TRUE, method = "jitter", jitter = 0.2)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

boxplot(SES$SES_fbsor ~ trmt, data = SES, col = c("gray60", "aquamarine3"),
        ylab = "Standardized Effect Sizes (SES)", xlab = "", main = "Functional", cex.main = 2,
        cex.lab = 1.6, cex.axis = 1.2, ylim = c(-12,5), outline = FALSE)
stripchart(SES$SES_fbsor ~ trmt, data = SES, col = c("gray28", "aquamarine4"), vertical = TRUE,
           pch = 19, cex = 2, add = TRUE, method = "jitter", jitter = 0.2)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

boxplot(SES$SES_pbsor ~ trmt, data = SES, col = c("gray60", "aquamarine3"),
        ylab = "Standardized Effect Sizes (SES)", xlab = "", main = "Phylogenetic", cex.main = 2,
        cex.lab = 1.6, cex.axis = 1.2, ylim = c(-12,5), outline = FALSE)
stripchart(SES$SES_pbsor ~ trmt, data = SES, col = c("gray28", "aquamarine4"), vertical = TRUE,
           pch = 19, cex = 2, add = TRUE, method = "jitter", jitter = 0.2)
abline(h = 0.0, col = "black", lwd = 3, lty=2)


dev.off()
