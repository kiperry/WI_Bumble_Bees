###################################################################################
#
# Wisconsin bumble bee data
#
# Step 2: Urban/Agricultural Pools to Local Greenspaces
#
# Community Trait Means
# Taxonomic Diversity
# Functional Diversity
# Phylogenetic Diversity
# Null Models
# PLS Analyses
# Landscape threshold analyses
#
# KI Perry; 21 December 2021
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
  spBB <- randomizeMatrix(samp = a, null.model = "independentswap")
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

write.csv(nqbl_var, file = "Urban_to_Local_2Axes/nqbl_var.u.csv")
write.csv(nmbl_var, file = "Urban_to_Local_2Axes/nmbl_var.u.csv")
write.csv(nwbl_var, file = "Urban_to_Local_2Axes/nwbl_var.u.csv")
write.csv(ntl_0, file = "Urban_to_Local_2Axes/ntl_0.u.csv")
write.csv(ntl_1, file = "Urban_to_Local_2Axes/ntl_1.u.csv")
write.csv(ntl_2, file = "Urban_to_Local_2Axes/ntl_2.u.csv")
write.csv(nnestl_0, file = "Urban_to_Local_2Axes/nnestl_0.u.csv")
write.csv(nnestl_2, file = "Urban_to_Local_2Axes/nnestl_2.u.csv")
write.csv(nit, file = "Urban_to_Local_2Axes/nit.u.csv")
write.csv(nwingl, file = "Urban_to_Local_2Axes/nwingl.u.csv")
write.csv(nheadw, file = "Urban_to_Local_2Axes/nheadw.u.csv")
write.csv(neyel, file = "Urban_to_Local_2Axes/neyel.u.csv")
write.csv(nthairl, file = "Urban_to_Local_2Axes/nthairl.u.csv")
write.csv(nbsetael, file = "Urban_to_Local_2Axes/nbsetael.u.csv")
write.csv(nbtl, file = "Urban_to_Local_2Axes/nbtl.u.csv")

write.csv(nbsim, file = "Urban_to_Local_2Axes/bb_tbeta_sim.u.csv")
write.csv(nbsne, file = "Urban_to_Local_2Axes/bb_tbeta_sne.u.csv")
write.csv(nbsor, file = "Urban_to_Local_2Axes/bb_tbeta_sor.u.csv")

write.csv(nfalpha, file = "Urban_to_Local_2Axes/bb_falpha.u.csv")
write.csv(nfsim, file = "Urban_to_Local_2Axes/bb_fbeta_sim.u.csv")
write.csv(nfsne, file = "Urban_to_Local_2Axes/bb_fbeta_sne.u.csv")
write.csv(nfsor, file = "Urban_to_Local_2Axes/bb_fbeta_sor.u.csv")

write.csv(npsim, file = "Urban_to_Local_2Axes/bb_pbeta_sim.csv")
write.csv(npsne, file = "Urban_to_Local_2Axes/bb_pbeta_sne.csv")
write.csv(npsor, file = "Urban_to_Local_2Axes/bb_pbeta_sor.csv")

############################################################################
# load the datasets
nqbl_var <- read.csv("Urban_to_Local_2Axes/nqbl_var.u.csv", row.names=1)
nmbl_var <- read.csv("Urban_to_Local_2Axes/nmbl_var.u.csv", row.names=1)
nwbl_var <- read.csv("Urban_to_Local_2Axes/nwbl_var.u.csv", row.names=1)
ntl_0 <- read.csv("Urban_to_Local_2Axes/ntl_0.u.csv", row.names=1)
ntl_1 <- read.csv("Urban_to_Local_2Axes/ntl_1.u.csv", row.names=1)
ntl_2 <- read.csv("Urban_to_Local_2Axes/ntl_2.u.csv", row.names=1)
nnestl_0 <- read.csv("Urban_to_Local_2Axes/nnestl_0.u.csv", row.names=1)
nnestl_2 <- read.csv("Urban_to_Local_2Axes/nnestl_2.u.csv", row.names=1)
nit <- read.csv("Urban_to_Local_2Axes/nit.u.csv", row.names=1)
nwingl <- read.csv("Urban_to_Local_2Axes/nwingl.u.csv", row.names=1)
nheadw <- read.csv("Urban_to_Local_2Axes/nheadw.u.csv", row.names=1)
neyel <- read.csv("Urban_to_Local_2Axes/neyel.u.csv", row.names=1)
nthairl <- read.csv("Urban_to_Local_2Axes/nthairl.u.csv", row.names=1)
nbsetael <- read.csv("Urban_to_Local_2Axes/nbsetael.u.csv", row.names=1)
nbtl <- read.csv("Urban_to_Local_2Axes/nbtl.u.csv", row.names=1)

nbsim <- read.csv("Urban_to_Local_2Axes/bb_tbeta_sim.u.csv", row.names=1)
nbsne <- read.csv("Urban_to_Local_2Axes/bb_tbeta_sne.u.csv", row.names=1)
nbsor <- read.csv("Urban_to_Local_2Axes/bb_tbeta_sor.u.csv", row.names=1)

nfalpha <- read.csv("Urban_to_Local_2Axes/bb_falpha.u.csv", row.names=1)
nfsim <- read.csv("Urban_to_Local_2Axes/bb_fbeta_sim.u.csv", row.names=1)
nfsne <- read.csv("Urban_to_Local_2Axes/bb_fbeta_sne.u.csv", row.names=1)
nfsor <- read.csv("Urban_to_Local_2Axes/bb_fbeta_sor.u.csv", row.names=1)

npsim<- read.csv(file = "Urban_to_Local_2Axes/bb_pbeta_sim.csv", row.names=1)
npsne<- read.csv(file = "Urban_to_Local_2Axes/bb_pbeta_sne.csv", row.names=1)
npsor<- read.csv(file = "Urban_to_Local_2Axes/bb_pbeta_sor.csv", row.names=1)

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

write.csv(SES, file = "SES_Local.csv")

# load the dataset
SES <- read.csv("SES_Local.csv", row.names=1)


###################################################################################
####Differences in Overall Observed SES from Null Communities######################


## pull out data for each treatment
ag <- SES[which(SES$trmt == "Agriculture"),]
str(ag)

urban <- SES[which(SES$trmt == "Urban"),]
str(urban)


## now let's compare by treatment

## Queen body length variance
hist(SES_qblv)
plot(SES_qblv, pch = 19, cex = 1.5)
abline(h = 0.0, col = "black", lwd = 3, lty=2)
boxplot(SES_qblv ~ SES$trmt)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

## compare to null expectations by treatment
w.qbl_var.u <- wilcox.test(urban$SES_qblv, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)#
w.qbl_var.u
w.qbl_var.a <- wilcox.test(ag$SES_qblv, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)#
w.qbl_var.a


## Male body length variance
hist(SES_mblv)
plot(SES_mblv, pch = 19, cex = 1.5)
abline(h = 0.0, col = "black", lwd = 3, lty=2)
boxplot(SES_mblv ~ SES$trmt)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

## compare to null expectations by treatment
w.mbl_var.u <- wilcox.test(urban$SES_mblv, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)#
w.mbl_var.u
w.mbl_var.a <- wilcox.test(ag$SES_mblv, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)#
w.mbl_var.a


## Worker body length variance
hist(SES_wblv)
plot(SES_wblv, pch = 19, cex = 1.5)
abline(h = 0.0, col = "black", lwd = 3, lty=2)
boxplot(SES_wblv ~ SES$trmt)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

## compare to null expectations by treatment
w.wbl_var.u <- wilcox.test(urban$SES_wblv, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)#
w.wbl_var.u
w.wbl_var.a <- wilcox.test(ag$SES_wblv, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)#
w.wbl_var.a


## short tongue
hist(SES_tl_0)
plot(SES_tl_0, pch = 19, cex = 1.5)
abline(h = 0.0, col = "black", lwd = 3, lty=2)
boxplot(SES_tl_0 ~ SES$trmt)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

## compare to null expectations by treatment
w.tl_0.u <- wilcox.test(urban$SES_tl_0, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.tl_0.u
w.tl_0.a <- wilcox.test(ag$SES_tl_0, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.tl_0.a


## medium tongue
hist(SES_tl_1)
plot(SES_tl_1, pch = 19, cex = 1.5, ylim = c(-0.5, 1.8))
abline(h = 0.0, col = "black", lwd = 3, lty=2)
boxplot(SES_tl_1 ~ SES$trmt, ylim = c(-0.5, 1.8))
abline(h = 0.0, col = "black", lwd = 3, lty=2)

## compare to null expectations by treatment
w.tl_1.u <- wilcox.test(urban$SES_tl_1, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.tl_1.u
w.tl_1.a <- wilcox.test(ag$SES_tl_1, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.tl_1.a


## long tongue
hist(SES_tl_2)
plot(SES_tl_2, pch = 19, cex = 1.5)
abline(h = 0.0, col = "black", lwd = 3, lty=2)
boxplot(SES_tl_2 ~ SES$trmt)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

## compare to null expectations by treatment
w.tl_2.u <- wilcox.test(urban$SES_tl_2, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.tl_2.u
w.tl_2.a <- wilcox.test(ag$SES_tl_2, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.tl_2.a


## below ground nests
hist(SES_nest0)
plot(SES_nest0, pch = 19, cex = 1.5)
abline(h = 0.0, col = "black", lwd = 3, lty=2)
boxplot(SES_nest0 ~ SES$trmt)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

## compare to null expectations by treatment
w.nest0.u <- wilcox.test(urban$SES_nest0, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.nest0.u
w.nest0.a <- wilcox.test(ag$SES_nest0, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.nest0.a


## above ground nests
hist(SES_nest2)
plot(SES_nest2, pch = 19, cex = 1.5)
abline(h = 0.0, col = "black", lwd = 3, lty=2)
boxplot(SES_nest2 ~ SES$trmt)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

## compare to null expectations by treatment
w.nest2.u <- wilcox.test(urban$SES_nest2, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.nest2.u
w.nest2.a <- wilcox.test(ag$SES_nest2, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.nest2.a


## inter-tegular distance (body size)
hist(SES_it)
plot(SES_it, pch = 19, cex = 1.5, ylim = c(-0.5, 16))
abline(h = 0.0, col = "black", lwd = 3, lty=2)
boxplot(SES_it ~ SES$trmt, ylim = c(-0.5, 16))
abline(h = 0.0, col = "black", lwd = 3, lty=2)

## compare to null expectations by treatment
w.it.u <- wilcox.test(urban$SES_it, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.it.u
w.it.a <- wilcox.test(ag$SES_it, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.it.a


## wing length
hist(SES_wingl)
plot(SES_wingl, pch = 19, cex = 1.5)
abline(h = 0.0, col = "black", lwd = 3, lty=2)
boxplot(SES_wingl ~ SES$trmt)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

## compare to null expectations by treatment
w.wingl.u <- wilcox.test(urban$SES_wingl, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.wingl.u
w.wingl.a <- wilcox.test(ag$SES_wingl, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.wingl.a


## head width
hist(SES_headw)
plot(SES_headw, pch = 19, cex = 1.5)
abline(h = 0.0, col = "black", lwd = 3, lty=2)
boxplot(SES_headw ~ SES$trmt)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

## compare to null expectations by treatment
w.headw.u <- wilcox.test(urban$SES_headw, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.headw.u
w.headw.a <- wilcox.test(ag$SES_headw, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.headw.a


## eye length
hist(SES_eyel)
plot(SES_eyel, pch = 19, cex = 1.5)
abline(h = 0.0, col = "black", lwd = 3, lty=2)
boxplot(SES_eyel ~ SES$trmt)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

## compare to null expectations by treatment
w.eyel.u <- wilcox.test(urban$SES_eyel, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.eyel.u
w.eyel.a <- wilcox.test(ag$SES_eyel, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.eyel.a


## thorax hair length
hist(SES_thairl)
plot(SES_thairl, pch = 19, cex = 1.5, ylim = c(-0.5, 45))
abline(h = 0.0, col = "black", lwd = 3, lty=2)
boxplot(SES_thairl ~ SES$trmt, ylim = c(-0.5, 45))
abline(h = 0.0, col = "black", lwd = 3, lty=2)

## compare to null expectations by treatment
w.thairl.u <- wilcox.test(urban$SES_thairl, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.thairl.u
w.thairl.a <- wilcox.test(ag$SES_thairl, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.thairl.a


## tibia setae length
hist(SES_setael)
plot(SES_setael, pch = 19, cex = 1.5)
abline(h = 0.0, col = "black", lwd = 3, lty=2)
boxplot(SES_setael ~ SES$trmt)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

## compare to null expectations by treatment
w.setael.u <- wilcox.test(urban$SES_setael, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.setael.u
w.setael.a <- wilcox.test(ag$SES_setael, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.setael.a


## tibia length
hist(SES_tibial)
plot(SES_tibial, pch = 19, cex = 1.5)
abline(h = 0.0, col = "black", lwd = 3, lty=2)
boxplot(SES_tibial ~ SES$trmt)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

## compare to null expectations by treatment
w.tibial.u <- wilcox.test(urban$SES_tibial, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.tibial.u
w.tibial.a <- wilcox.test(ag$SES_tibial, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.tibial.a


## taxonomic diveristy - beta sor
hist(SES_bsor)
plot(SES_bsor, pch = 19, cex = 1.5)
abline(h = 0.0, col = "black", lwd = 3, lty=2)
boxplot(SES_bsor ~ SES$trmt)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

## compare to null expectations by treatment
w.bsor.u <- wilcox.test(urban$SES_bsor, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.bsor.u
w.bsor.a <- wilcox.test(ag$SES_bsor, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.bsor.a


## taxonomic diveristy - beta sim
hist(SES_bsim)
plot(SES_bsim, pch = 19, cex = 1.5)
abline(h = 0.0, col = "black", lwd = 3, lty=2)
boxplot(SES_bsim ~ SES$trmt)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

## compare to null expectations by treatment
w.bsim.u <- wilcox.test(urban$SES_bsim, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.bsim.u
w.bsim.a <- wilcox.test(ag$SES_bsim, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.bsim.a


## taxonomic diveristy - beta sne
hist(SES_bsne)
plot(SES_bsne, pch = 19, cex = 1.5)
abline(h = 0.0, col = "black", lwd = 3, lty=2)
boxplot(SES_bsne ~ SES$trmt)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

## compare to null expectations by treatment
w.bsne.u <- wilcox.test(urban$SES_bsne, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.bsne.u
w.bsne.a <- wilcox.test(ag$SES_bsne, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.bsne.a


## functional alpha diversity
hist(SES_falpha)
plot(SES_falpha, pch = 19, cex = 1.5)
abline(h = 0.0, col = "black", lwd = 3, lty=2)
boxplot(SES_falpha ~ SES$trmt)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

## compare to null expectations by treatment
w.falpha.u <- wilcox.test(urban$SES_falpha, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.falpha.u
w.falpha.a <- wilcox.test(ag$SES_falpha, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.falpha.a


## functional beta diversity - beta sor
hist(SES_fbsor)
plot(SES_fbsor, pch = 19, cex = 1.5)
abline(h = 0.0, col = "black", lwd = 3, lty=2)
boxplot(SES_fbsor ~ SES$trmt)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

## compare to null expectations by treatment
w.fbsor.u <- wilcox.test(urban$SES_fbsor, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.fbsor.u
w.fbsor.a <- wilcox.test(ag$SES_fbsor, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.fbsor.a


## functional beta diversity - beta sim
hist(SES_fbsim)
plot(SES_fbsim, pch = 19, cex = 1.5)
abline(h = 0.0, col = "black", lwd = 3, lty=2)
boxplot(SES_fbsim ~ SES$trmt)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

## compare to null expectations by treatment
w.fbsim.u <- wilcox.test(urban$SES_fbsim, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.fbsim.u
w.fbsim.a <- wilcox.test(ag$SES_fbsim, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.fbsim.a


## functional beta diversity - beta sne
hist(SES_fbsne)
plot(SES_fbsne, pch = 19, cex = 1.5)
abline(h = 0.0, col = "black", lwd = 3, lty=2)
boxplot(SES_fbsne ~ SES$trmt)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

## compare to null expectations by treatment
w.fbsne.u <- wilcox.test(urban$SES_fbsne, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.fbsne.u
w.fbsne.a <- wilcox.test(ag$SES_fbsne, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.fbsne.a


## phylogenetic beta diversity - beta sor
hist(SES_pbsor)
plot(SES_pbsor, pch = 19, cex = 1.5)
abline(h = 0.0, col = "black", lwd = 3, lty=2)
boxplot(SES_pbsor ~ SES$trmt)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

## compare to null expectations by treatment
w.pbsor.u <- wilcox.test(urban$SES_pbsor, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.pbsor.u
w.pbsor.a <- wilcox.test(ag$SES_pbsor, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.pbsor.a


## phylogenetic beta diversity - beta sim
hist(SES_pbsim)
plot(SES_pbsim, pch = 19, cex = 1.5)
abline(h = 0.0, col = "black", lwd = 3, lty=2)
boxplot(SES_pbsim ~ SES$trmt)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

## compare to null expectations by treatment
w.pbsim.u <- wilcox.test(urban$SES_pbsim, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.pbsim.u
w.pbsim.a <- wilcox.test(ag$SES_pbsim, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.pbsim.a


## phylogenetic beta diversity - beta sne
hist(SES_pbsne)
plot(SES_pbsne, pch = 19, cex = 1.5)
abline(h = 0.0, col = "black", lwd = 3, lty=2)
boxplot(SES_pbsne ~ SES$trmt)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

## compare to null expectations by treatment
w.pbsne.u <- wilcox.test(urban$SES_pbsne, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.pbsne.u
w.pbsne.a <- wilcox.test(ag$SES_pbsne, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.pbsne.a


########################################################################################################
## Figures

# Taxonomic beta-diversity
png("L_SES_TDiv.png", width = 2500, height = 1000, pointsize = 30)

par(mfrow=c(1,3))
par(mar=c(5,8,4,2))

boxplot(SES$SES_bsor ~ trmt, data = SES, col = c("gray60", "aquamarine3"),
        ylab = "Standardized Effect Sizes (SES)", xlab = "", main = "(a) Total Beta-Diversity", cex.main = 2,
        cex.lab = 1.6, cex.axis = 1.2, ylim = c(-2.5,2.5), outline = FALSE)
stripchart(SES$SES_bsor ~ trmt, data = SES, col = c("gray28", "aquamarine4"), vertical = TRUE,
           pch = 19, cex = 2, add = TRUE, method = "jitter", jitter = 0.2)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

boxplot(SES$SES_bsim ~ trmt, data = SES, col = c("gray60", "aquamarine3"),
        ylab = "Standardized Effect Sizes (SES)", xlab = "", main = "(b) Turnover", cex.main = 2,
        cex.lab = 1.6, cex.axis = 1.2, ylim = c(-2.5,2.5), outline = FALSE)
stripchart(SES$SES_bsim ~ trmt, data = SES, col = c("gray28", "aquamarine4"), vertical = TRUE,
           pch = 19, cex = 2, add = TRUE, method = "jitter", jitter = 0.2)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

boxplot(SES$SES_bsne ~ trmt, data = SES, col = c("gray60", "aquamarine3"),
        ylab = "Standardized Effect Sizes (SES)", xlab = "", main = "(c) Nestedness", cex.main = 2,
        cex.lab = 1.6, cex.axis = 1.2, ylim = c(-2.5,2.5), outline = FALSE)
stripchart(SES$SES_bsne ~ trmt, data = SES, col = c("gray28", "aquamarine4"), vertical = TRUE,
           pch = 19, cex = 2, add = TRUE, method = "jitter", jitter = 0.2)
abline(h = 0.0, col = "black", lwd = 3, lty=2)


dev.off()


# Functional beta-diversity
png("L_SES_FDiv.png", width = 2500, height = 1000, pointsize = 30)

par(mfrow=c(1,3))
par(mar=c(5,8,4,2))

boxplot(SES$SES_fbsor ~ trmt, data = SES, col = c("gray60", "aquamarine3"),
        ylab = "Standardized Effect Sizes (SES)", xlab = "", main = "(a) Total Beta-Diversity", cex.main = 2,
        cex.lab = 1.6, cex.axis = 1.2, ylim = c(-1,2.5), outline = FALSE)
stripchart(SES$SES_fbsor ~ trmt, data = SES, col = c("gray28", "aquamarine4"), vertical = TRUE,
           pch = 19, cex = 2, add = TRUE, method = "jitter", jitter = 0.2)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

boxplot(SES$SES_fbsim ~ trmt, data = SES, col = c("gray60", "aquamarine3"),
        ylab = "Standardized Effect Sizes (SES)", xlab = "", main = "(b) Turnover", cex.main = 2,
        cex.lab = 1.6, cex.axis = 1.2, ylim = c(-1,2.5), outline = FALSE)
stripchart(SES$SES_fbsim ~ trmt, data = SES, col = c("gray28", "aquamarine4"), vertical = TRUE,
           pch = 19, cex = 2, add = TRUE, method = "jitter", jitter = 0.2)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

boxplot(SES$SES_fbsne ~ trmt, data = SES, col = c("gray60", "aquamarine3"),
        ylab = "Standardized Effect Sizes (SES)", xlab = "", main = "(c) Nestedness", cex.main = 2,
        cex.lab = 1.6, cex.axis = 1.2, ylim = c(-1,2.5), outline = FALSE)
stripchart(SES$SES_fbsne ~ trmt, data = SES, col = c("gray28", "aquamarine4"), vertical = TRUE,
           pch = 19, cex = 2, add = TRUE, method = "jitter", jitter = 0.2)
abline(h = 0.0, col = "black", lwd = 3, lty=2)


dev.off()


# Phylogenetic beta-diversity
png("L_SES_PDiv.png", width = 2500, height = 1000, pointsize = 30)

par(mfrow=c(1,3))
par(mar=c(5,8,4,2))

boxplot(SES$SES_pbsor ~ trmt, data = SES, col = c("gray60", "aquamarine3"),
        ylab = "Standardized Effect Sizes (SES)", xlab = "", main = "(a) Total Beta-Diversity", cex.main = 2,
        cex.lab = 1.6, cex.axis = 1.2, ylim = c(-1,4), outline = FALSE)
stripchart(SES$SES_pbsor ~ trmt, data = SES, col = c("gray28", "aquamarine4"), vertical = TRUE,
           pch = 19, cex = 2, add = TRUE, method = "jitter", jitter = 0.2)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

boxplot(SES$SES_pbsim ~ trmt, data = SES, col = c("gray60", "aquamarine3"),
        ylab = "Standardized Effect Sizes (SES)", xlab = "", main = "(b) Turnover", cex.main = 2,
        cex.lab = 1.6, cex.axis = 1.2, ylim = c(-1,4), outline = FALSE)
stripchart(SES$SES_pbsim ~ trmt, data = SES, col = c("gray28", "aquamarine4"), vertical = TRUE,
           pch = 19, cex = 2, add = TRUE, method = "jitter", jitter = 0.2)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

boxplot(SES$SES_pbsne ~ trmt, data = SES, col = c("gray60", "aquamarine3"),
        ylab = "Standardized Effect Sizes (SES)", xlab = "", main = "(c) Nestedness", cex.main = 2,
        cex.lab = 1.6, cex.axis = 1.2, ylim = c(-1,4), outline = FALSE)
stripchart(SES$SES_pbsne ~ trmt, data = SES, col = c("gray28", "aquamarine4"), vertical = TRUE,
           pch = 19, cex = 2, add = TRUE, method = "jitter", jitter = 0.2)
abline(h = 0.0, col = "black", lwd = 3, lty=2)


dev.off()

#####################################################################################################
### figure panel - all diversity indices

library(ggplot2)
library(ggthemes)
library(reshape2)
library(viridis)

SES_tbeta <- as.data.frame(SES[,c(18,17,16)])
colnames(SES_tbeta) <- c("Nestedness", "Turnover", "Total")
SES_tbeta$trmt <- c("U", "A", "U", "U", "U", "A", "A", "U", "A", "U", "U", "U", "U", "A", "U",
                    "A", "U", "U", "A", "A")
SES_tbeta$trmt <- as.factor(SES_tbeta$trmt)
str(SES_tbeta)

SES_fbeta <- as.data.frame(SES[,c(21,20,19)])
colnames(SES_fbeta) <- c("Nestedness", "Turnover", "Total")
SES_fbeta$trmt <- c("U", "A", "U", "U", "U", "A", "A", "U", "A", "U", "U", "U", "U", "A", "U",
                    "A", "U", "U", "A", "A")
SES_fbeta$trmt <- as.factor(SES_fbeta$trmt)
str(SES_fbeta)

SES_pbeta <- as.data.frame(SES[,c(25,24,23)])
colnames(SES_pbeta) <- c("Nestedness", "Turnover", "Total")
SES_pbeta$trmt <- c("U", "A", "U", "U", "U", "A", "A", "U", "A", "U", "U", "U", "U", "A", "U",
                    "A", "U", "U", "A", "A")
SES_pbeta$trmt <- as.factor(SES_pbeta$trmt)
str(SES_pbeta)

SES_tbeta$metric <- rep(c("Taxonomic"),each = 20)
SES_fbeta$metric <- rep(c("Functional"),each = 20)
SES_pbeta$metric <- rep(c("Phylogenetic"),each = 20)

SES_div <- rbind(SES_tbeta, SES_fbeta, SES_pbeta)
SES_div.m <- melt(SES_div)
colnames(SES_div.m) <- c("trmt","metric", "var", "ses")

levels(SES_div.m$trmt)[levels(SES_div.m$trmt)=="U"] <- "Urban"
levels(SES_div.m$trmt)[levels(SES_div.m$trmt)=="A"] <- "Agriculture"

SES_div.m$metric <- factor(SES_div.m$metric, levels = c("Taxonomic", "Functional", "Phylogenetic"))


png("SES_Div_Local.png", width = 2000, height = 1000, pointsize = 20)

ggplot(SES_div.m, aes(x=ses, y=var, fill = var)) +
  geom_boxplot(outlier.shape = NA) +
  geom_dotplot(position = position_jitter(width = 0.2, height = 0.2),
               dotsize = 2,
               binaxis = "y",
               stackdir = "center") +
  facet_grid(metric ~ trmt) + 
  coord_cartesian(xlim = c(-1, 3)) +
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


##############################################################################

# load the dataset
SES <- read.csv("SES_Local.csv", row.names=1)

## import the environmental data

## landscape

land <- read.csv("./landscape_bytransect.csv", row.names=1)
str(land)

## check for correlated variables
## pull out variables to include in the analysis
plot(land[2:8], pch = 19)
cor(land[2:8], method = c("pearson"), use = "complete.obs")


## check variables for outliers

dotchart(land$SIDI1500)
plot(land$SIDI1500)

dotchart(land$pland.ag.1500)
plot(land$pland.ag.1500)

dotchart(land$pland.ur.1500)
plot(land$pland.ur.1500)

dotchart(land$pland.nat.1500)
plot(land$pland.nat.1500)

dotchart(land$lpi.for.1500)
plot(land$lpi.for.1500)

dotchart(land$ed.for.1500)
plot(land$ed.for.1500)

dotchart(land$enn.for.1500)
plot(land$enn.for.1500)

##########################################################################################
##PLS analysis for pooled (2019 and 2020) bumble bee data and landscape data
##SES values

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("mixOmics")
library(mixOmics)

## landscape

## community trait means
cwm_land <- pls(land[2:8], SES[1:15], mode = c("canonical"), ncomp = 2, scale = TRUE, max.iter = 100)
cwm_land

cwm_land$loadings.star
cwm_land$prop_expl_var
cwm_land$loadings #use 0.3 as cutoff
cwm_land$variates

plotVar(cwm_land)
plotLoadings(cwm_land)

## run a reduced pls with variables that made the initial 0.3 cutoff
cwm_land.red <- pls(land[2:7], SES[,c(2,4,6,9,11:13,15)], mode = c("canonical"), ncomp = 2, scale = TRUE, max.iter = 100)
cwm_land.red
cwm_land.red$loadings
cwm_land.red$loadings.star
cwm_land.red$prop_expl_var

plotVar(cwm_land.red)
plotLoadings(cwm_land.red)
cim(cwm_land.red)$mat.cor

nw_cwm_land <- network(cwm_land.red, cutoff = 0.5, color.edge = color.spectral(2), lty.edge = c("solid", "dashed"),
            lwd.edge = 2)
nw_cwm_land

cwm_nw <- as.data.frame(nw_cwm_land$M)
write.csv(cwm_nw, file = "cwm_similarity_values.csv")

cwm_land.red.pred <- as.data.frame(cwm_land.red$loadings$X)
write.csv(cwm_land.red.pred, file = "cwm_predictor_loadings.csv")

cwm_land.red.resp <- as.data.frame(cwm_land.red$loadings$Y)
write.csv(cwm_land.red.resp, file = "cwm_response_loadings.csv")

plot(cwm_land.red.pred$comp2 ~ cwm_land.red.pred$comp1, pch = 19, ylim = c(-0.7, 0.7), xlim = c(-0.7, 0.7), col = "gray47",
     xlab = "PLS Axis 1", ylab = "PLS Axis 2", main = "(a) CWM Metrics", cex = 1.6)
points(cwm_land.red.resp$comp2 ~ cwm_land.red.resp$comp1, pch = 15, cex = 1.6)
abline(h = 0.0, v = 0.0, col = "black", lwd = 1, lty=1)
text(-0.52, 0.402, "Landscape Diversity", pos = 4, font = 1, cex = 1)
text(0.04, -0.63, "Agriculture", pos = 4, font = 1, cex = 1)
text(0.19, 0.621, "Urban", pos = 4, font = 1, cex = 1)
text(-0.48, -0.185, "Natural Habitat", pos = 4, font = 1, cex = 1)
text(-0.42, 0.069, "LPI Forest", pos = 4, font = 1, cex = 1)
text(-0.53, -0.10, "ED Forest", pos = 2, font = 1, cex = 1)

text(-0.43, 0.185, "MBLV", pos = 2, font = 2, cex = 1)
text(0.05, 0.647, "Short Tongue", pos = 2, font = 2, cex = 1)
text(-0.148, -0.577, "Long Tongue", pos = 2, font = 2, cex = 1)
text(-0.33, -0.08, "IT", pos = 2, font = 2, cex = 1)
text(0.49, 0.064, "Head Width", pos = 2, font = 2, cex = 1)
text(0.45, 0.197, "Eye Length", pos = 2, font = 2, cex = 1)
text(0.36, -0.188, "Hair Length", pos = 2, font = 2, cex = 1)
text(0.31, -0.375, "Tibia Length", pos = 2, font = 2, cex = 1)

## diversity metrics
div_land <- pls(land[2:8], SES[16:25], mode = c("canonical"), ncomp = 2, scale = TRUE, max.iter = 100)
div_land

div_land$prop_expl_var
div_land$loadings #use 0.3 as cutoff 
div_land$loadings.star
div_land$variates

plotVar(div_land)
plotLoadings(div_land)

## run a reduced pls with variables that made the initial 0.3 cutoff
div_land.red <- pls(land[2:8], SES[,c(16:19,22:25)], mode = c("canonical"), ncomp = 2, scale = TRUE, max.iter = 100)
div_land.red
div_land$prop_expl_var
div_land$loadings
div_land$loadings.star

plotVar(div_land.red)
plotLoadings(div_land.red)
cim(div_land.red)$mat.cor

nw_div_land <- network(div_land.red, cutoff = 0.5, color.edge = color.spectral(2), lty.edge = c("solid", "dashed"),
                       lwd.edge = 2)
nw_div_land

div_nw <- as.data.frame(nw_div_land$M)
write.csv(div_nw, file = "div_similarity_values.csv")

div_land.red.pred <- as.data.frame(div_land.red$loadings$X)
write.csv(cwm_land.red.pred, file = "div_predictor_loadings.csv")

div_land.red.resp <- as.data.frame(div_land.red$loadings$Y)
write.csv(cwm_land.red.resp, file = "div_response_loadings.csv")

plot(div_land.red.pred$comp2 ~ div_land.red.pred$comp1, pch = 19, ylim = c(-0.7, 0.7), xlim = c(-0.7, 0.7), col = "gray47",
     xlab = "PLS Axis 1", ylab = "PLS Axis 2", main = "(a) Diversity Metrics", cex = 1.6)
points(div_land.red.resp$comp2 ~ div_land.red.resp$comp1, pch = 15, cex = 1.6)
abline(h = 0.0, v = 0.0, col = "black", lwd = 1, lty=1)
text(0.69, -0.225, "Landscape Diversity", pos = 2, font = 1, cex = 1)
text(-0.52, -0.278, "Agriculture", pos = 4, font = 1, cex = 1)
text(0.45, 0.296, "Urban", pos = 4, font = 1, cex = 1)
text(0.04, -0.36, "Natural Habitat", pos = 4, font = 1, cex = 1)
text(0.17, -0.30, "LPI Forest", pos = 3, font = 1, cex = 1)
text(0.029, -0.54, "ED Forest", pos = 4, font = 1, cex = 1)
text(-0.09, 0.518, "ENN Forest", pos = 2, font = 1, cex = 1)

text(-0.395, -0.128, "TBsor", pos = 2, font = 2, cex = 1)
text(-0.38, -0.115, "TBsim", pos = 4, font = 2, cex = 1)
text(0.40, 0.069, "TBsne", pos = 2, font = 2, cex = 1)
text(-0.065, 0.363, "FBsor", pos = 2, font = 2, cex = 1)
text(0.35, -0.31, "Falpha", pos = 2, font = 2, cex = 1)
text(-0.347, -0.568, "PBsor", pos = 2, font = 2, cex = 1)
text(-0.39, -0.11, "PBsim", pos = 3, font = 2, cex = 1)
text(0.375, -0.62, "PBsne", pos = 2, font = 2, cex = 1)


# Make the figure
png("Fig.PLS.png", width = 2200, height = 1000, pointsize = 20)

par(mfrow=c(1,2)) # one row and two columns
par(mar=c(5,8,4,2))

#CWM
plot(cwm_land.red.pred$comp2 ~ cwm_land.red.pred$comp1, pch = 19, ylim = c(-0.7, 0.7), xlim = c(-0.7, 0.7), col = "gray47",
     xlab = "PLS Axis 1", ylab = "PLS Axis 2", main = "(a) CWM Metrics", cex = 1.6)
points(cwm_land.red.resp$comp2 ~ cwm_land.red.resp$comp1, pch = 15, cex = 1.6)
abline(h = 0.0, v = 0.0, col = "black", lwd = 1, lty=1)
text(-0.52, 0.402, "Landscape Diversity", pos = 4, font = 1, cex = 1)
text(0.04, -0.63, "Agriculture", pos = 4, font = 1, cex = 1)
text(0.19, 0.621, "Urban", pos = 4, font = 1, cex = 1)
text(-0.48, -0.185, "Natural Habitat", pos = 4, font = 1, cex = 1)
text(-0.42, 0.069, "LPI Forest", pos = 4, font = 1, cex = 1)
text(-0.53, -0.10, "ED Forest", pos = 2, font = 1, cex = 1)

text(-0.43, 0.185, "MBLV", pos = 2, font = 2, cex = 1)
text(0.05, 0.644, "Short Tongue", pos = 2, font = 2, cex = 1)
text(-0.148, -0.577, "Long Tongue", pos = 2, font = 2, cex = 1)
text(-0.33, -0.08, "IT", pos = 2, font = 2, cex = 1)
text(0.49, 0.064, "Head Width", pos = 2, font = 2, cex = 1)
text(0.45, 0.197, "Eye Length", pos = 2, font = 2, cex = 1)
text(0.36, -0.188, "Hair Length", pos = 2, font = 2, cex = 1)
text(0.31, -0.375, "Corbicula Length", pos = 2, font = 2, cex = 1)

#Diversity
plot(div_land.red.pred$comp2 ~ div_land.red.pred$comp1, pch = 19, ylim = c(-0.7, 0.7), xlim = c(-0.7, 0.7), col = "gray47",
     xlab = "PLS Axis 1", ylab = "PLS Axis 2", main = "(b) Diversity Metrics", cex = 1.6)
points(div_land.red.resp$comp2 ~ div_land.red.resp$comp1, pch = 15, cex = 1.6)
abline(h = 0.0, v = 0.0, col = "black", lwd = 1, lty=1)
text(0.69, -0.225, "Landscape Diversity", pos = 2, font = 1, cex = 1)
text(-0.52, -0.278, "Agriculture", pos = 4, font = 1, cex = 1)
text(0.45, 0.296, "Urban", pos = 4, font = 1, cex = 1)
text(0.04, -0.36, "Natural Habitat", pos = 4, font = 1, cex = 1)
text(0.17, -0.30, "LPI Forest", pos = 3, font = 1, cex = 1)
text(0.029, -0.54, "ED Forest", pos = 4, font = 1, cex = 1)
text(-0.09, 0.518, "ENN Forest", pos = 2, font = 1, cex = 1)

text(-0.395, -0.128, "TBsor", pos = 2, font = 2, cex = 1)
text(-0.38, -0.115, "TBsim", pos = 4, font = 2, cex = 1)
text(0.40, 0.069, "TBsne", pos = 2, font = 2, cex = 1)
text(-0.065, 0.363, "FBsor", pos = 2, font = 2, cex = 1)
text(0.35, -0.31, "Falpha", pos = 2, font = 2, cex = 1)
text(-0.347, -0.568, "PBsor", pos = 2, font = 2, cex = 1)
text(-0.39, -0.11, "PBsim", pos = 3, font = 2, cex = 1)
text(0.375, -0.627, "PBsne", pos = 2, font = 2, cex = 1)

dev.off()
