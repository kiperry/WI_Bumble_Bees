###################################################################################
#
# Wisconsin bumble bee data
#
# Step 2: Urban/Agricultural Pools to Local Greenspaces
# Revised analysis that removed the trait shuffle randomization from the null model
#
# Community Trait Means
# Taxonomic Diversity
# Functional Diversity
# Phylogenetic Diversity
# Null Models
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

############################################################################
# load the datasets
nqbl_var <- read.csv("Urban_to_Local_2Axes_NoTraitShuffle/nqbl_var.u.csv", row.names=1)
nmbl_var <- read.csv("Urban_to_Local_2Axes_NoTraitShuffle/nmbl_var.u.csv", row.names=1)
nwbl_var <- read.csv("Urban_to_Local_2Axes_NoTraitShuffle/nwbl_var.u.csv", row.names=1)
ntl_0 <- read.csv("Urban_to_Local_2Axes_NoTraitShuffle/ntl_0.u.csv", row.names=1)
ntl_1 <- read.csv("Urban_to_Local_2Axes_NoTraitShuffle/ntl_1.u.csv", row.names=1)
ntl_2 <- read.csv("Urban_to_Local_2Axes_NoTraitShuffle/ntl_2.u.csv", row.names=1)
nnestl_0 <- read.csv("Urban_to_Local_2Axes_NoTraitShuffle/nnestl_0.u.csv", row.names=1)
nnestl_2 <- read.csv("Urban_to_Local_2Axes_NoTraitShuffle/nnestl_2.u.csv", row.names=1)
nit <- read.csv("Urban_to_Local_2Axes_NoTraitShuffle/nit.u.csv", row.names=1)
nwingl <- read.csv("Urban_to_Local_2Axes_NoTraitShuffle/nwingl.u.csv", row.names=1)
nheadw <- read.csv("Urban_to_Local_2Axes_NoTraitShuffle/nheadw.u.csv", row.names=1)
neyel <- read.csv("Urban_to_Local_2Axes_NoTraitShuffle/neyel.u.csv", row.names=1)
nthairl <- read.csv("Urban_to_Local_2Axes_NoTraitShuffle/nthairl.u.csv", row.names=1)
nbsetael <- read.csv("Urban_to_Local_2Axes_NoTraitShuffle/nbsetael.u.csv", row.names=1)
nbtl <- read.csv("Urban_to_Local_2Axes_NoTraitShuffle/nbtl.u.csv", row.names=1)

nbsim <- read.csv("Urban_to_Local_2Axes_NoTraitShuffle/bb_tbeta_sim.u.csv", row.names=1)
nbsne <- read.csv("Urban_to_Local_2Axes_NoTraitShuffle/bb_tbeta_sne.u.csv", row.names=1)
nbsor <- read.csv("Urban_to_Local_2Axes_NoTraitShuffle/bb_tbeta_sor.u.csv", row.names=1)

nfalpha <- read.csv("Urban_to_Local_2Axes_NoTraitShuffle/bb_falpha.u.csv", row.names=1)
nfsim <- read.csv("Urban_to_Local_2Axes_NoTraitShuffle/bb_fbeta_sim.u.csv", row.names=1)
nfsne <- read.csv("Urban_to_Local_2Axes_NoTraitShuffle/bb_fbeta_sne.u.csv", row.names=1)
nfsor <- read.csv("Urban_to_Local_2Axes_NoTraitShuffle/bb_fbeta_sor.u.csv", row.names=1)

npsim<- read.csv(file = "Urban_to_Local_2Axes_NoTraitShuffle/bb_pbeta_sim.csv", row.names=1)
npsne<- read.csv(file = "Urban_to_Local_2Axes_NoTraitShuffle/bb_pbeta_sne.csv", row.names=1)
npsor<- read.csv(file = "Urban_to_Local_2Axes_NoTraitShuffle/bb_pbeta_sor.csv", row.names=1)


#################################################################################################################

## SES Calculations
#calculate standardized effect sizes (SES) for each trait and index
#the effect size is the difference between the observed value and the expected one
#then divide the effect size by the standard deviation of the null distribution to get the standardized effect size
#allows comparison among sites with different numbers of species


##calculate SES values for each metric

## community weighted means
## queen body length variance
SES_qblv <- (cwm.obs$qbl_var - apply(nqbl_var, MARGIN = 1, mean)) / apply(nqbl_var, MARGIN = 1, sd, na.rm=T)
SES_qblv

## male body length variance
SES_mblv <- (cwm.obs$mbl_var - apply(nmbl_var, MARGIN = 1, mean)) / apply(nmbl_var, MARGIN = 1, sd, na.rm=T)
SES_mblv

## worker body length variance
SES_wblv <- (cwm.obs$wbl_var - apply(nwbl_var, MARGIN = 1, mean)) / apply(nwbl_var, MARGIN = 1, sd, na.rm=T)
SES_wblv

## tongue length - short
SES_tl_0 <- (cwm.obs$tl_0 - apply(ntl_0, MARGIN = 1, mean)) / apply(ntl_0, MARGIN = 1, sd, na.rm=T)
SES_tl_0

## tongue length - medium
SES_tl_1 <- (cwm.obs$tl_1 - apply(ntl_1, MARGIN = 1, mean)) / apply(ntl_1, MARGIN = 1, sd, na.rm=T)
SES_tl_1

## tongue length - long
SES_tl_2 <- (cwm.obs$tl_2 - apply(ntl_2, MARGIN = 1, mean)) / apply(ntl_2, MARGIN = 1, sd, na.rm=T)
SES_tl_2

## nest location - belowground
SES_nest0 <- (cwm.obs$nestl_0 - apply(nnestl_0, MARGIN = 1, mean)) / apply(nnestl_0, MARGIN = 1, sd, na.rm=T)
SES_nest0

## nest location - aboveground
SES_nest2 <- (cwm.obs$nestl_2 - apply(nnestl_2, MARGIN = 1, mean)) / apply(nnestl_2, MARGIN = 1, sd, na.rm=T)
SES_nest2

## intertegular distance
SES_it <- (cwm.obs$it - apply(nit, MARGIN = 1, mean)) / apply(nit, MARGIN = 1, sd, na.rm=T)
SES_it

## wing length
SES_wingl <- (cwm.obs$rwingl - apply(nwingl, MARGIN = 1, mean)) / apply(nwingl, MARGIN = 1, sd, na.rm=T)
SES_wingl

## head width
SES_headw <- (cwm.obs$rheadw - apply(nheadw, MARGIN = 1, mean)) / apply(nheadw, MARGIN = 1, sd, na.rm=T)
SES_headw

## eye length
SES_eyel <- (cwm.obs$reyel - apply(neyel, MARGIN = 1, mean)) / apply(neyel, MARGIN = 1, sd, na.rm=T)
SES_eyel

## thorax hair length
SES_thairl <- (cwm.obs$thairl - apply(nthairl, MARGIN = 1, mean)) / apply(nthairl, MARGIN = 1, sd, na.rm=T)
SES_thairl

## corbicula setae length
SES_setael <- (cwm.obs$bsetael - apply(nbsetael, MARGIN = 1, mean)) / apply(nbsetael, MARGIN = 1, sd, na.rm=T)
SES_setael

## corbicula length
SES_tibial <- (cwm.obs$rbtl - apply(nbtl, MARGIN = 1, mean)) / apply(nbtl, MARGIN = 1, sd, na.rm=T)
SES_tibial


## Diversity Indices##

#Taxonomic diversity

beta.sor <- as.matrix(bb.dist$beta.sor)
beta.sor <- colMeans(beta.sor)

beta.sim <- as.matrix(bb.dist$beta.sim)
beta.sim <- colMeans(beta.sim)

beta.sne <- as.matrix(bb.dist$beta.sne)
beta.sne <- colMeans(beta.sne)

beta.t <- data.frame(beta.sor, beta.sim, beta.sne)

## taxonomic diveristy - beta sor
SES_bsor <- (beta.t$beta.sor - apply(nbsor, MARGIN = 1, mean)) / apply(nbsor, MARGIN = 1, sd, na.rm=T)
SES_bsor

## taxonomic diveristy - beta sim
SES_bsim <- (beta.t$beta.sim - apply(nbsim, MARGIN = 1, mean)) / apply(nbsim, MARGIN = 1, sd, na.rm=T)
SES_bsim

## taxonomic diveristy - beta sne
SES_bsne <- (beta.t$beta.sne - apply(nbsne, MARGIN = 1, mean)) / apply(nbsne, MARGIN = 1, sd, na.rm=T)
SES_bsne


# Functional Diversity
falpha <- as.matrix(bb.falpha)

fbeta.sor <- as.matrix(bb.fun$funct.beta.sor)
fbeta.sor <- colMeans(fbeta.sor)

fbeta.sim <- as.matrix(bb.fun$funct.beta.sim)
fbeta.sim <- colMeans(fbeta.sim)

fbeta.sne <- as.matrix(bb.fun$funct.beta.sne)
fbeta.sne <- colMeans(fbeta.sne)

beta.f <- data.frame(fbeta.sor, fbeta.sim, fbeta.sne, falpha)

## functional alpha diversity
SES_falpha <- (beta.f$falpha - apply(nfalpha, MARGIN = 1, mean)) / apply(nfalpha, MARGIN = 1, sd, na.rm=T)
SES_falpha

## functional diveristy - beta sor
SES_fbsor <- (beta.f$fbeta.sor - apply(nfsor, MARGIN = 1, mean)) / apply(nfsor, MARGIN = 1, sd, na.rm=T)
SES_fbsor

## functional diveristy - beta sim
SES_fbsim <- (beta.f$fbeta.sim - apply(nfsim, MARGIN = 1, mean)) / apply(nfsim, MARGIN = 1, sd, na.rm=T)
SES_fbsim

## functional diveristy - beta sne
SES_fbsne <- (beta.f$fbeta.sne - apply(nfsne, MARGIN = 1, mean)) / apply(nfsne, MARGIN = 1, sd, na.rm=T)
SES_fbsne


# Phylogenetic Diversity
pbeta.sor <- as.matrix(bb.phy$phylo.beta.sor)
pbeta.sor <- colMeans(pbeta.sor)

pbeta.sim <- as.matrix(bb.phy$phylo.beta.sim)
pbeta.sim <- colMeans(pbeta.sim)

pbeta.sne <- as.matrix(bb.phy$phylo.beta.sne)
pbeta.sne <- colMeans(pbeta.sne)

beta.p <- data.frame(pbeta.sor, pbeta.sim, pbeta.sne)


## phylogenetic diversity - beta sor
SES_pbsor <- (beta.p$pbeta.sor - apply(npsor, MARGIN = 1, mean)) / apply(npsor, MARGIN = 1, sd, na.rm=T)
SES_pbsor

## phylogenetic diversity - beta sim
SES_pbsim <- (beta.p$pbeta.sim - apply(npsim, MARGIN = 1, mean)) / apply(npsim, MARGIN = 1, sd, na.rm=T)
SES_pbsim

## phylogenetic diversity - beta sne
SES_pbsne <- (beta.p$pbeta.sne - apply(npsne, MARGIN = 1, mean)) / apply(npsne, MARGIN = 1, sd, na.rm=T)
SES_pbsne


## combine all indices into one matrix
SES <- cbind(SES_qblv, SES_mblv, SES_wblv, SES_tl_0, SES_tl_1, SES_tl_2, SES_nest0,
             SES_nest2, SES_it, SES_wingl, SES_headw, SES_eyel, SES_thairl, SES_setael, SES_tibial,
             SES_bsor, SES_bsim, SES_bsne, SES_fbsor, SES_fbsim, SES_fbsne, SES_falpha, SES_pbsor, 
             SES_pbsim, SES_pbsne)
SES <- as.data.frame(SES)
str(SES)

SES$trmt <- c("Urban", "Agriculture", "Urban", "Urban", "Urban", "Agriculture", "Agriculture",
              "Urban", "Agriculture", "Urban", "Urban", "Urban", "Urban", "Agriculture", "Urban",
              "Agriculture", "Urban", "Urban", "Agriculture", "Agriculture")
str(SES)
SES$trmt <- as.factor(SES$trmt)
str(SES)

write.csv(SES, file = "SES_Local_NoTraitShuffleModel.csv")

# load the dataset
SES <- read.csv("SES_Local_NoTraitShuffleModel.csv", row.names=1)


###################################################################################
####Differences in Overall Observed SES from Null Communities######################


## pull out data for each treatment
ag <- SES[which(SES$trmt == "Agriculture"),]
str(ag)

urban <- SES[which(SES$trmt == "Urban"),]
str(urban)


## now let's compare by treatment

## Queen body length variance
hist(SES$SES_qblv)
plot(SES$SES_qblv, pch = 19, cex = 1.5)
abline(h = 0.0, col = "black", lwd = 3, lty=2)
boxplot(SES$SES_qblv ~ SES$trmt)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

## compare to null expectations by treatment
w.qbl_var.u <- wilcox.test(urban$SES_qblv, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)#
w.qbl_var.u
w.qbl_var.a <- wilcox.test(ag$SES_qblv, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)#
w.qbl_var.a


## Male body length variance
hist(SES$SES_mblv)
plot(SES$SES_mblv, pch = 19, cex = 1.5)
abline(h = 0.0, col = "black", lwd = 3, lty=2)
boxplot(SES$SES_mblv ~ SES$trmt)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

## compare to null expectations by treatment
w.mbl_var.u <- wilcox.test(urban$SES_mblv, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)#
w.mbl_var.u
w.mbl_var.a <- wilcox.test(ag$SES_mblv, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)#
w.mbl_var.a


## Worker body length variance
hist(SES$SES_wblv)
plot(SES$SES_wblv, pch = 19, cex = 1.5)
abline(h = 0.0, col = "black", lwd = 3, lty=2)
boxplot(SES$SES_wblv ~ SES$trmt)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

## compare to null expectations by treatment
w.wbl_var.u <- wilcox.test(urban$SES_wblv, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)#
w.wbl_var.u
w.wbl_var.a <- wilcox.test(ag$SES_wblv, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)#
w.wbl_var.a


## short tongue
hist(SES$SES_tl_0)
plot(SES$SES_tl_0, pch = 19, cex = 1.5)
abline(h = 0.0, col = "black", lwd = 3, lty=2)
boxplot(SES$SES_tl_0 ~ SES$trmt)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

## compare to null expectations by treatment
w.tl_0.u <- wilcox.test(urban$SES_tl_0, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.tl_0.u
w.tl_0.a <- wilcox.test(ag$SES_tl_0, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.tl_0.a


## medium tongue
hist(SES$SES_tl_1)
plot(SES$SES_tl_1, pch = 19, cex = 1.5, ylim = c(-0.5, 1.8))
abline(h = 0.0, col = "black", lwd = 3, lty=2)
boxplot(SES$SES_tl_1 ~ SES$trmt, ylim = c(-0.5, 1.8))
abline(h = 0.0, col = "black", lwd = 3, lty=2)

## compare to null expectations by treatment
w.tl_1.u <- wilcox.test(urban$SES_tl_1, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.tl_1.u
w.tl_1.a <- wilcox.test(ag$SES_tl_1, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.tl_1.a


## long tongue
hist(SES$SES_tl_2)
plot(SES$SES_tl_2, pch = 19, cex = 1.5)
abline(h = 0.0, col = "black", lwd = 3, lty=2)
boxplot(SES$SES_tl_2 ~ SES$trmt)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

## compare to null expectations by treatment
w.tl_2.u <- wilcox.test(urban$SES_tl_2, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.tl_2.u
w.tl_2.a <- wilcox.test(ag$SES_tl_2, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.tl_2.a


## below ground nests
hist(SES$SES_nest0)
plot(SES$SES_nest0, pch = 19, cex = 1.5)
abline(h = 0.0, col = "black", lwd = 3, lty=2)
boxplot(SES$SES_nest0 ~ SES$trmt)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

## compare to null expectations by treatment
w.nest0.u <- wilcox.test(urban$SES_nest0, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.nest0.u
w.nest0.a <- wilcox.test(ag$SES_nest0, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.nest0.a


## above ground nests
hist(SES$SES_nest2)
plot(SES$SES_nest2, pch = 19, cex = 1.5)
abline(h = 0.0, col = "black", lwd = 3, lty=2)
boxplot(SES$SES_nest2 ~ SES$trmt)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

## compare to null expectations by treatment
w.nest2.u <- wilcox.test(urban$SES_nest2, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.nest2.u
w.nest2.a <- wilcox.test(ag$SES_nest2, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.nest2.a


## inter-tegular distance (body size)
hist(SES$SES_it)
plot(SES$SES_it, pch = 19, cex = 1.5)
abline(h = 0.0, col = "black", lwd = 3, lty=2)
boxplot(SES$SES_it ~ SES$trmt)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

## compare to null expectations by treatment
w.it.u <- wilcox.test(urban$SES_it, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.it.u
w.it.a <- wilcox.test(ag$SES_it, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.it.a


## wing length
hist(SES$SES_wingl)
plot(SES$SES_wingl, pch = 19, cex = 1.5)
abline(h = 0.0, col = "black", lwd = 3, lty=2)
boxplot(SES$SES_wingl ~ SES$trmt)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

## compare to null expectations by treatment
w.wingl.u <- wilcox.test(urban$SES_wingl, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.wingl.u
w.wingl.a <- wilcox.test(ag$SES_wingl, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.wingl.a


## head width
hist(SES$SES_headw)
plot(SES$SES_headw, pch = 19, cex = 1.5)
abline(h = 0.0, col = "black", lwd = 3, lty=2)
boxplot(SES$SES_headw ~ SES$trmt)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

## compare to null expectations by treatment
w.headw.u <- wilcox.test(urban$SES_headw, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.headw.u
w.headw.a <- wilcox.test(ag$SES_headw, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.headw.a


## eye length
hist(SES$SES_eyel)
plot(SES$SES_eyel, pch = 19, cex = 1.5)
abline(h = 0.0, col = "black", lwd = 3, lty=2)
boxplot(SES$SES_eyel ~ SES$trmt)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

## compare to null expectations by treatment
w.eyel.u <- wilcox.test(urban$SES_eyel, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.eyel.u
w.eyel.a <- wilcox.test(ag$SES_eyel, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.eyel.a


## thorax hair length
hist(SES$SES_thairl)
plot(SES$SES_thairl, pch = 19, cex = 1.5)
abline(h = 0.0, col = "black", lwd = 3, lty=2)
boxplot(SES$SES_thairl ~ SES$trmt)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

## compare to null expectations by treatment
w.thairl.u <- wilcox.test(urban$SES_thairl, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.thairl.u
w.thairl.a <- wilcox.test(ag$SES_thairl, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.thairl.a


## tibia setae length
hist(SES$SES_setael)
plot(SES$SES_setael, pch = 19, cex = 1.5)
abline(h = 0.0, col = "black", lwd = 3, lty=2)
boxplot(SES$SES_setael ~ SES$trmt)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

## compare to null expectations by treatment
w.setael.u <- wilcox.test(urban$SES_setael, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.setael.u
w.setael.a <- wilcox.test(ag$SES_setael, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.setael.a


## tibia length
hist(SES$SES_tibial)
plot(SES$SES_tibial, pch = 19, cex = 1.5)
abline(h = 0.0, col = "black", lwd = 3, lty=2)
boxplot(SES$SES_tibial ~ SES$trmt)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

## compare to null expectations by treatment
w.tibial.u <- wilcox.test(urban$SES_tibial, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.tibial.u
w.tibial.a <- wilcox.test(ag$SES_tibial, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.tibial.a


## taxonomic diveristy - beta sor
hist(SES$SES_bsor)
plot(SES$SES_bsor, pch = 19, cex = 1.5)
abline(h = 0.0, col = "black", lwd = 3, lty=2)
boxplot(SES$SES_bsor ~ SES$trmt)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

## compare to null expectations by treatment
w.bsor.u <- wilcox.test(urban$SES_bsor, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.bsor.u
w.bsor.a <- wilcox.test(ag$SES_bsor, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.bsor.a


## taxonomic diveristy - beta sim
hist(SES$SES_bsim)
plot(SES$SES_bsim, pch = 19, cex = 1.5)
abline(h = 0.0, col = "black", lwd = 3, lty=2)
boxplot(SES$SES_bsim ~ SES$trmt)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

## compare to null expectations by treatment
w.bsim.u <- wilcox.test(urban$SES_bsim, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.bsim.u
w.bsim.a <- wilcox.test(ag$SES_bsim, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.bsim.a


## taxonomic diveristy - beta sne
hist(SES$SES_bsne)
plot(SES$SES_bsne, pch = 19, cex = 1.5)
abline(h = 0.0, col = "black", lwd = 3, lty=2)
boxplot(SES$SES_bsne ~ SES$trmt)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

## compare to null expectations by treatment
w.bsne.u <- wilcox.test(urban$SES_bsne, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.bsne.u
w.bsne.a <- wilcox.test(ag$SES_bsne, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.bsne.a


## functional alpha diversity
hist(SES$SES_falpha)
plot(SES$SES_falpha, pch = 19, cex = 1.5)
abline(h = 0.0, col = "black", lwd = 3, lty=2)
boxplot(SES$SES_falpha ~ SES$trmt)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

## compare to null expectations by treatment
w.falpha.u <- wilcox.test(urban$SES_falpha, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.falpha.u
w.falpha.a <- wilcox.test(ag$SES_falpha, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.falpha.a


## functional beta diversity - beta sor
hist(SES$SES_fbsor)
plot(SES$SES_fbsor, pch = 19, cex = 1.5)
abline(h = 0.0, col = "black", lwd = 3, lty=2)
boxplot(SES$SES_fbsor ~ SES$trmt)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

## compare to null expectations by treatment
w.fbsor.u <- wilcox.test(urban$SES_fbsor, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.fbsor.u
w.fbsor.a <- wilcox.test(ag$SES_fbsor, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.fbsor.a


## functional beta diversity - beta sim
hist(SES$SES_fbsim)
plot(SES$SES_fbsim, pch = 19, cex = 1.5)
abline(h = 0.0, col = "black", lwd = 3, lty=2)
boxplot(SES$SES_fbsim ~ SES$trmt)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

## compare to null expectations by treatment
w.fbsim.u <- wilcox.test(urban$SES_fbsim, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.fbsim.u
w.fbsim.a <- wilcox.test(ag$SES_fbsim, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.fbsim.a


## functional beta diversity - beta sne
hist(SES$SES_fbsne)
plot(SES$SES_fbsne, pch = 19, cex = 1.5)
abline(h = 0.0, col = "black", lwd = 3, lty=2)
boxplot(SES$SES_fbsne ~ SES$trmt)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

## compare to null expectations by treatment
w.fbsne.u <- wilcox.test(urban$SES_fbsne, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.fbsne.u
w.fbsne.a <- wilcox.test(ag$SES_fbsne, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.fbsne.a


## phylogenetic beta diversity - beta sor
hist(SES$SES_pbsor)
plot(SES$SES_pbsor, pch = 19, cex = 1.5)
abline(h = 0.0, col = "black", lwd = 3, lty=2)
boxplot(SES$SES_pbsor ~ SES$trmt)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

## compare to null expectations by treatment
w.pbsor.u <- wilcox.test(urban$SES_pbsor, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.pbsor.u
w.pbsor.a <- wilcox.test(ag$SES_pbsor, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.pbsor.a


## phylogenetic beta diversity - beta sim
hist(SES$SES_pbsim)
plot(SES$SES_pbsim, pch = 19, cex = 1.5)
abline(h = 0.0, col = "black", lwd = 3, lty=2)
boxplot(SES$SES_pbsim ~ SES$trmt)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

## compare to null expectations by treatment
w.pbsim.u <- wilcox.test(urban$SES_pbsim, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.pbsim.u
w.pbsim.a <- wilcox.test(ag$SES_pbsim, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.pbsim.a


## phylogenetic beta diversity - beta sne
hist(SES$SES_pbsne)
plot(SES$SES_pbsne, pch = 19, cex = 1.5)
abline(h = 0.0, col = "black", lwd = 3, lty=2)
boxplot(SES$SES_pbsne ~ SES$trmt)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

## compare to null expectations by treatment
w.pbsne.u <- wilcox.test(urban$SES_pbsne, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.pbsne.u
w.pbsne.a <- wilcox.test(ag$SES_pbsne, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.pbsne.a
