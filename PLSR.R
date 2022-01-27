###################################################################################
#
# Wisconsin bumble bee data
#
# PLSR analyses
#
# KI Perry; 6 August 2021
#
###################################################################################

t <- read.csv("./bb_rtraits.csv", row.names=1)
a <- read.csv("./bb_abund.csv", row.names=1)

str(a)
rowSums(a)
colSums(a)

# create a vector with the column sums for each species
# species not found in Cleveland will have a 0
sp <- colSums(a)
sp

# removes any columns (i.e. species) that are not found in Cleveland
a <- a[, colSums(a != 0) > 0]
str(a)
colSums(a)

#Convert abundance data to relative abundance
library(vegan)

#Check total abundance in each sample (i.e. row)
apply(a, 1, sum)

#Convert to relative abundance
ra <- decostand(a, method = "total")

#Check the totals for each row in the new dataset
#All rows should equal 1
apply(ra, 1, sum)

str(ra)

# add the sp vector as a column in the trait matrix, shows which species
# are found in Cleveland and which are absent (i.e. with a 0)
t$sp <- sp
t

# uses the sp values to remove rows of species not collected in Cleveland
# then remove the column because we don't need it anymore
t <- t[t$sp != 0, ]
t <- t[,-28]
str(t)

#trim the trait dataset by identifying traits that are highly correlated
#or lack sufficient variance among species
names(t)
str(t)
plot(t)
cor(t, method = c("pearson"), use = "complete.obs")

#lecty, nest construction, sociality, activity, parasitism, and pollen transport lack sufficient
#variance among bumble bee species - remove these traits
t2 <- t[,-8]#lecty
t2 <- t2[,-9]#nest construction
t2 <- t2[,-9]#sociality
t2 <- t2[,-9]#activity
t2 <- t2[,-9]#parasitism
t2 <- t2[,-9]#pollen transport

t2 <- t2[,-21]#wt

plot(t2)
cor(t2, method = c("pearson"), use = "complete.obs")

t2 <- t2[,-20]#corbicula width

t2 <- t2[,-10]#wing marginal cell length (keeping inter-tegular distance)
t2 <- t2[,-10]#wing width

t2 <- t2[,-12]#eye width

plot(t2)
cor(t2, method = c("pearson"), use = "complete.obs")

#remove the body size variables from the literature except body length variances
#these are correlated with each other and several other traits
t2 <- t2[,-1]#average queen body length
t2 <- t2[,-2]#average male body length
t2 <- t2[,-3]#average worker body length
t2 <- t2[,-10] #scape length
t2 <- t2[,-12] #tibia length
t2 <- t2[,-8] #head width

plot(t2, pch = 19)
cor(t2, method = c("pearson"), use = "complete.obs")

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
t3$wbl_var <- log(t3$wbl_var + 1)

hist(t2$it)
hist(log(t2$it))#
t3$it <- log(t3$it + 1)

hist(t2$rwingl)
hist(log(t2$rwingl))
t3$rwingl <- log(t3$rwingl + 1)

hist(t2$reyel)
hist(log(t2$reyel) + 1)#
t3$reyel <- log(t3$reyel + 1)

hist(t2$thairl)
hist(log(t2$thairl))#
t3$thairl <- log(t3$thairl + 1)

hist(t2$bsetael)
hist(log(t2$bsetael))#


#Double check all species present in both datasets
intersect(colnames(a), rownames(t3))
names(a)
#Double check if a species is present in one dataset but not the other
setdiff(colnames(a), rownames(t3))
setdiff(rownames (t3), colnames(a))

rownames(t3) == colnames(a) # we are good!


##############################################################################
## Community Metrics

library(FD)
library(picante)
library(gawdis)
library(betapart)

#observed CWM
cwm.obs <- functcomp(t3, as.matrix(ra), CWM.type = "all")
cwm.obs

#observed taxonomic beta diversity
# create beta part object for analyses
str(ra)
bb.core <- betapart.core.abund(ra)

# returns three dissimilarity matrices containing 
# pairwise between-site values of each beta-diversity component
bb.dist <- beta.pair.abund(bb.core, index.family = "bray")
str(bb.dist)

#observed functional beta diversity
plot(t3, pch = 19)
str(t3)
cor(t3[6:10], method = c("pearson"), use = "complete.obs")

tdis <- gawdis(t3, w.type = "optimized", opti.maxiter = 300,
               groups.weight = T, groups = c(1, 2, 3, 4, 5, 6, 6, 6, 6, 6))
attr(tdis, "correls")
attr(tdis, "weights")

#group head traits that are correlated to down weight them in distance matrix
pcoBB <- dudi.pco(sqrt(tdis), scannf = FALSE, nf = 4)#select four axes
scatter(pcoBB)

pcoBB$li
sum(pcoBB$eig[1:3]) / sum(pcoBB$eig)
sum(pcoBB$eig[1:2]) / sum(pcoBB$eig)

# check correlations among axes and traits
str(t2)
cor(pcoBB$li, t2, use = "complete.obs")

#due to number of bees at each site, can only use first two axes of PCoA
t.ax <- as.matrix(pcoBB$li[1:2])
#bb.fcore <- functional.betapart.core(a, t4, warning.time = FALSE)takes too long to run
bb.fun <- functional.beta.pair(ra, t.ax, index.family = "sorensen")
str(bb.fun)

#observed functional alpha diversity
bb.rao <- Rao(sample = t(a), dfunc = tdis, dphyl = NULL, weight = FALSE, Jost = TRUE, structure = NULL)
bb.falpha <- bb.rao$FD$Alpha
bb.falpha


##############################################################################

## import the environmental data

land.all <- read.csv("./landscape_bytransect.csv", row.names=1)
str(land.all)
land <- land.all[1:17] #focus only on 1500 m
str(land)

## check for correlated variables
## pull out variables to include in the analysis
plot(land[2:17], pch = 19)
cor(land[2:17], method = c("pearson"), use = "complete.obs")


land <- as.data.frame(cbind(land$SIDI1500, land$pland.ag.1500, land$pland.ur.1500, land$pland.nat.1500,
                            land$lpi.for.1500, land$clumpy.for.1500, land$enn.for.1500))

colnames(land) <- c("SIDI1500", "pland.ag.1500", "pland.ur.1500", "pland.nat.1500", "lpi.for.1500",
                    "clumpy.for.1500", "enn.for.1500")
plot(land, pch = 19)
cor(land, method = c("pearson"), use = "complete.obs")

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

dotchart(land$clumpy.for.1500)
plot(land$clumpy.for.1500)
land <- land[,-6]
str(land)

dotchart(land$lpi.for.1500)
plot(land$lpi.for.1500)


##########################################################################################
##PLSR analysis for pooled (2019 and 2020) bumble bee data
#install.packages("plsdepot")
library(plsdepot)

## cwm
cwm_1500 <- plsreg2(land, cwm.obs, comps = 2, crosval = TRUE)
plot(cwm_1500)

cwm_1500$cor.xt##these are correct for variables
cwm_1500$cor.yu##these are correct for taxa

cwm_1500$std.coefs
cwm_1500$resid
cwm_1500$expvar
cwm_1500$VIP
cwm_1500$Q2cum


cwm_500 <- plsreg2(land500, SES[1:15], comps = 2, crosval = TRUE)
plot(cwm_500)

cwm_500$cor.xt##these are correct for variables
cwm_500$cor.yu##these are correct for taxa

cwm_500$std.coefs
cwm_500$resid
cwm_500$expvar
cwm_500$VIP
cwm_500$Q2cum


## diversity indices
div_1500 <- plsreg2(land15, SES[16:22], comps = 2, crosval = TRUE)
plot(div_1500)

div_1500$cor.xt##these are correct for variables
div_1500$cor.yu##these are correct for taxa

div_1500$std.coefs
div_1500$resid
div_1500$expvar
div_1500$VIP
div_1500$Q2cum

div_500 <- plsreg2(land500, SES[16:22], comps = 2, crosval = TRUE)
plot(div_500)

div_500$cor.xt##these are correct for variables
div_500$cor.yu##these are correct for taxa

div_500$std.coefs
div_500$resid
div_500$expvar
div_500$VIP
div_500$Q2cum

