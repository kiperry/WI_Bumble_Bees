###################################################################################
#
# Wisconsin bumble bee data
#
# RLQ and fourth corner analyses
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

## import the environmental data

land.all <- read.csv("./landscape_bytransect.csv", row.names=1)
str(land.all)
land <- land[1:17] #focus only on 1500 m
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


library(ade4)

## RLQ - landscape
##RLQ analysis for pooled (2019 and 2020) bumble bee data
abund.L.land <- dudi.coa(a, scannf = FALSE)
env.R.land <- dudi.pca(land, row.w = abund.L.land$lw, scannf = FALSE, center = TRUE, scale = TRUE)
trait.Q.land <- dudi.hillsmith(t3, row.w = abund.L.land$cw, scannf = FALSE)

rlq.bb.land <- rlq(env.R.land, abund.L.land, trait.Q.land, scannf = FALSE, nf = 3)

summary(rlq.bb.land)
print(rlq.bb.land)

plot(rlq.bb.land)
par(mfrow = c(1, 3))
s.arrow(rlq.bb.land$l1)
s.arrow(rlq.bb.land$c1)
s.label(rlq.bb.land$lQ, boxes = FALSE)

rlq.bb.land$eig[1]/sum(rlq.bb.land$eig)
rlq.bb.land$eig[2]/sum(rlq.bb.land$eig)

## weighted correlations axes / env.
t(env.R.land$tab)%*%(diag(env.R.land$lw))%*%as.matrix(rlq.bb.land$mR)

## weighted correlations axes / traits.
t(trait.Q.land$tab)%*%(diag(trait.Q.land$lw))%*%as.matrix(rlq.bb.land$mQ)

##A biplot representing traits and environmental variables
par(mfrow = c(1, 1))
s.arrow(rlq.bb.land$c1, xlim=c(-1.5,1.5), boxes = FALSE)
s.label(rlq.bb.land$l1, add.plot=T, clab=1.0)

## correlations traits / env.
rlq.bb.land$tab
rlq.bb.land$aR #Corr env.R.land axes / rlq axes
rlq.bb.land$aQ #Corr abund.L.land axes / coinertia axes
rlq.bb.land$l1 #Principal components (loadings for env.R.land cols)
rlq.bb.land$c1 #Principal axes (loadings for trait.Q.land cols) 
rlq.bb.land$li #CT row scores (cols of env.R.land) 
rlq.bb.land$co #CT col scores (cols of trait.Q.land)


## based on principal component loadings - trim trait and landscape datasets
## use 0.4 cutoff

land <- land[,-2]
land <- land[,-4]
land <- land[,-4]
str(land)

t4 <- as.data.frame(cbind(t3$mbl_var, t3$wbl_var, t3$tl, t3$nestl, t3$it))
colnames(t4) <- c("mbl_var", "wbl_var", "tl", "nestl", "it")
t4$tl <- as.factor(t4$tl)
t4$nestl <- as.factor(t4$nestl)
str(t4)

t4 <- t4[,-4]


## RLQ with reduced datasets
abund.L.land <- dudi.coa(a, scannf = FALSE)
env.R.land <- dudi.pca(land, row.w = abund.L.land$lw, scannf = FALSE, center = TRUE, scale = TRUE)
trait.Q.land <- dudi.hillsmith(t4, row.w = abund.L.land$cw, scannf = FALSE)

rlq.bb.land <- rlq(env.R.land, abund.L.land, trait.Q.land, scannf = FALSE, nf = 3)

summary(rlq.bb.land)
print(rlq.bb.land)

plot(rlq.bb.land)
par(mfrow = c(1, 3))
s.arrow(rlq.bb.land$l1)
s.arrow(rlq.bb.land$c1)
s.label(rlq.bb.land$lQ, boxes = FALSE)

rlq.bb.land$eig[1]/sum(rlq.bb.land$eig)
rlq.bb.land$eig[2]/sum(rlq.bb.land$eig)

## weighted correlations axes / env.
t(env.R.land$tab)%*%(diag(env.R.land$lw))%*%as.matrix(rlq.bb.land$mR)

## weighted correlations axes / traits.
t(trait.Q.land$tab)%*%(diag(trait.Q.land$lw))%*%as.matrix(rlq.bb.land$mQ)

##A biplot representing traits and environmental variables
par(mfrow = c(1, 1))
s.arrow(rlq.bb.land$c1, xlim=c(-1.5,1.5), boxes = FALSE)
s.label(rlq.bb.land$l1, add.plot=T, clab=1.0)

## correlations traits / env.
rlq.bb.land$tab
rlq.bb.land$aR #Corr env.R.land axes / rlq axes
rlq.bb.land$aQ #Corr abund.L.land axes / coinertia axes
rlq.bb.land$l1 #Principal components (loadings for env.R.land cols)
rlq.bb.land$c1 #Principal axes (loadings for trait.Q.land cols) 
rlq.bb.land$li #CT row scores (cols of env.R.land) 
rlq.bb.land$co #CT col scores (cols of trait.Q.land)


## Fourth-corner analysis - landscape
nrepet <- 999

##P-values not adjusted for multiple comparisons
four.comb.bb.land <- fourthcorner(land, a, t4, modeltype = 6, p.adjust.method.G = "none",
                                 p.adjust.method.D = "none", nrepet = nrepet)

plot(four.comb.bb.land, alpha = 0.1, stat = "D2")
four.comb.bb.land


##adjust P-values for multiple comparisons
four.comb.bb.land.adj <- p.adjust.4thcorner(four.comb.bb.land, p.adjust.method.G = "fdr", 
                                            p.adjust.method.D = "fdr")

plot(four.comb.bb.land.adj, alpha = 0.1, stat = "D2")

## Blue = negative relationships
## Red = positive relationships


##Combining both approaches
##Evaluates the global significance of the traits-environment relationships, based on total
##inertia of the RLQ analysis
testrlq.bb.land <- randtest(rlq.bb.land, modeltype = 6, nrepet = nrepet)
testrlq.bb.land
plot(testrlq.bb.land)

##Calculate the total inertia of RLQ analysis (Srlq)
Srlq.bb.land <- fourthcorner2(land, a, t3, modeltype = 6, p.adjust.method.G = "fdr", nrepet = nrepet)
Srlq.bb.land$trRLQ

plot(four.comb.bb.land.adj, x.rlq = rlq.bb.land, alpha = 0.1, stat = "D2", type = "biplot")


##tests the links between RLQ axes and traits (Qaxes) or environmental variables (Raxes)
testQaxes.bb.land <- fourthcorner.rlq(rlq.bb.land, modeltype = 6, typetest = "Q.axes", nrepet = nrepet,
                                      p.adjust.method.G = "fdr", p.adjust.method.D = "fdr")
print(testQaxes.bb.land, stat = "D")


testRaxes.bb.land <- fourthcorner.rlq(rlq.bb.land, modeltype = 6, typetest = "R.axes", nrepet = nrepet,
                                      p.adjust.method.G = "fdr", p.adjust.method.D = "fdr")
print(testRaxes.bb.land, stat = "D")


##plot results
par(mfrow = c(1, 2))
plot(testQaxes.bb.land, alpha = 0.1, type = "table", stat = "D2")
plot(testRaxes.bb.land, alpha = 0.1, type = "table", stat = "D2")

par(mfrow = c(1, 2))
plot(testQaxes.bb.land, alpha = 0.05, type = "biplot",
     stat = "D2", col = c("black", "blue", "orange", "green"))
plot(testRaxes.bb.land, alpha = 0.05, type = "biplot",
     stat = "D2", col = c("black", "blue", "orange", "green"))



