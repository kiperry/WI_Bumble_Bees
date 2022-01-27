###################################################################################
#
# Wisconsin bumble bee data
#
# Exploratory graphing of measured continuous traits
#
# KI Perry; 29 June 2021
#
###################################################################################

#The formula for calculating standard error
se <- function(x) sd(x, na.rm=TRUE)/sqrt(sum(!is.na(x)))

#Explore measured traits by species
t2 <- read.csv("./bb_traits_sp.csv")

t2$species
table(t2$species)#Number of individuals of each species

library(viridis)

#####Box plots with data points added

#all traits correlated with inter-tegular distance except hair lengths
#so use relative measurements
cor(t2[,2:23], method = c("pearson"), use = "complete.obs")
plot(t2[,2:23])

##intertegular distance
par(mar = c(5,9,4,2))
boxplot(it ~ species, data = t2, col = rainbow(16, alpha = 0.3),
        xlab = "Intertegular Distance (mm)", ylab = "",
        horizontal = TRUE, las = 1, range = 0, main = "Inter-tegular Distance")
stripchart(it ~ species, data = t2, col = rainbow(16),
           pch = 19, cex = 2, las = 1, add = TRUE,
           method = "jitter", jitter = 0.2)

png("Inter-tegular_Distance.png", width = 2200, height = 1000, pointsize = 20)
par(mfrow=c(1,1))
par(mar = c(5,9,4,2))
boxplot(it ~ species, data = t2, col = viridis(16, alpha = 0.5),
        xlab = "Inter-tegular Distance (mm)", ylab = "",
        horizontal = TRUE, las = 1, range = 0, main = "",
        cex.main = 2, cex.lab = 1.5, cex.axis = 1.2)
stripchart(it ~ species, data = t2, col = viridis(16),
           pch = 19, cex = 2, las = 1, add = TRUE,
           method = "jitter", jitter = 0.2)
dev.off()

tapply(t2$it, list(t2$species), mean)
tapply(t2$it, list(t2$species), se)


##marginal cell length
par(mar = c(5,9,4,2))
boxplot(wingmc ~ species, data = t2, col = rainbow(16, alpha = 0.3),
        xlab = "Marginal Cell Length (mm)", ylab = "",
        horizontal = TRUE, las = 1, range = 0, main = "Marginal Cell Length")
stripchart(wingmc ~ species, data = t2, col = rainbow(16),
           pch = 19, cex = 2, las = 1, add = TRUE,
           method = "jitter", jitter = 0.2)

tapply(t2$wingmc, list(t2$species), mean)
tapply(t2$wingmc, list(t2$species), se)

boxplot(rwingmc ~ species, data = t2, col = rainbow(16, alpha = 0.3),
        xlab = "Relative Marginal Cell Length", ylab = "",
        horizontal = TRUE, las = 1, range = 0, main = "Marginal Cell Length")
stripchart(rwingmc ~ species, data = t2, col = rainbow(16),
           pch = 19, cex = 2, las = 1, add = TRUE,
           method = "jitter", jitter = 0.2)


##wing width
par(mar = c(5,9,4,2))
boxplot(rwingw ~ species, data = t2, col = rainbow(16, alpha = 0.3),
        xlab = "Relative Wing Width", ylab = "",
        horizontal = TRUE, las = 1, range = 0, main = "Wing Width")
stripchart(rwingw ~ species, data = t2, col = rainbow(16),
           pch = 19, cex = 2, las = 1, add = TRUE,
           method = "jitter", jitter = 0.2)


##wing length
par(mar = c(5,9,4,2))
boxplot(rwingl ~ species, data = t2, col = rainbow(16, alpha = 0.3),
        xlab = "Relative Wing Length", ylab = "",
        horizontal = TRUE, las = 1, range = 0, main = "Wing Length")
stripchart(rwingl ~ species, data = t2, col = rainbow(16),
           pch = 19, cex = 2, las = 1, add = TRUE,
           method = "jitter", jitter = 0.2)


##head width
par(mar = c(5,9,4,2))
boxplot(rheadw ~ species, data = t2, col = rainbow(16, alpha = 0.3),
        xlab = "Relative Head Width", ylab = "",
        horizontal = TRUE, las = 1, range = 0, main = "Head Width")
stripchart(rheadw ~ species, data = t2, col = rainbow(16),
           pch = 19, cex = 2, las = 1, add = TRUE,
           method = "jitter", jitter = 0.2)


##eye width
par(mar = c(5,9,4,2))
boxplot(reyew ~ species, data = t2, col = rainbow(16, alpha = 0.3),
        xlab = "Relative Eye Width", ylab = "",
        horizontal = TRUE, las = 1, range = 0, main = "Eye Width")
stripchart(reyew ~ species, data = t2, col = rainbow(16),
           pch = 19, cex = 2, las = 1, add = TRUE,
           method = "jitter", jitter = 0.2)


##eye length
par(mar = c(5,9,4,2))
boxplot(reyel ~ species, data = t2, col = rainbow(16, alpha = 0.3),
        xlab = "Relative Eye Length", ylab = "",
        horizontal = TRUE, las = 1, range = 0, main = "Eye Length")
stripchart(reyel ~ species, data = t2, col = rainbow(16),
           pch = 19, cex = 2, las = 1, add = TRUE,
           method = "jitter", jitter = 0.2)


##scape length
par(mar = c(5,9,4,2))
boxplot(rscapel ~ species, data = t2, col = rainbow(16, alpha = 0.3),
        xlab = "Relative Scape Length", ylab = "",
        horizontal = TRUE, las = 1, range = 0, main = "Scape Length")
stripchart(rscapel ~ species, data = t2, col = rainbow(16),
           pch = 19, cex = 2, las = 1, add = TRUE,
           method = "jitter", jitter = 0.2)


##thorax hair length
par(mar = c(5,9,4,2))
boxplot(thairl ~ species, data = t2, col = rainbow(16, alpha = 0.3),
        xlab = "Thorax Hair Length (mm)", ylab = "",
        horizontal = TRUE, las = 1, range = 0, main = "Thorax Hair Length")
stripchart(thairl ~ species, data = t2, col = rainbow(16),
           pch = 19, cex = 2, las = 1, add = TRUE,
           method = "jitter", jitter = 0.2)

png("Thorax_Hair_Length.png", width = 2200, height = 1000, pointsize = 20)
par(mfrow=c(1,1))
par(mar = c(5,9,4,2))
boxplot(thairl ~ species, data = t2, col = viridis(16, alpha = 0.5),
        xlab = "Thorax Hair Length (mm)", ylab = "",
        horizontal = TRUE, las = 1, range = 0, main = "",
        cex.main = 2, cex.lab = 1.5, cex.axis = 1.2)
stripchart(thairl ~ species, data = t2, col = viridis(16),
           pch = 19, cex = 2, las = 1, add = TRUE,
           method = "jitter", jitter = 0.2)
dev.off()


##corbicula hair length
par(mar = c(5,9,4,2))
boxplot(bsetael ~ species, data = t2, col = rainbow(16, alpha = 0.3),
        xlab = "Corbicula Setae Length (mm)", ylab = "",
        horizontal = TRUE, las = 1, range = 0, main = "Setae Length")
stripchart(bsetael ~ species, data = t2, col = rainbow(16),
           pch = 19, cex = 2, las = 1, add = TRUE,
           method = "jitter", jitter = 0.2)


##corbicula length
par(mar = c(5,9,4,2))
boxplot(rbtl ~ species, data = t2, col = rainbow(16, alpha = 0.3),
        xlab = "Relative Corbicula Length", ylab = "",
        horizontal = TRUE, las = 1, range = 0, main = "Corbicula Length")
stripchart(rbtl ~ species, data = t2, col = rainbow(16),
           pch = 19, cex = 2, las = 1, add = TRUE,
           method = "jitter", jitter = 0.2)


##corbicula width
par(mar = c(5,9,4,2))
boxplot(rbtw ~ species, data = t2, col = rainbow(16, alpha = 0.3),
        xlab = "Relative Corbicula Width", ylab = "",
        horizontal = TRUE, las = 1, range = 0, main = "Corbicula Width")
stripchart(rbtw ~ species, data = t2, col = rainbow(16),
           pch = 19, cex = 2, las = 1, add = TRUE,
           method = "jitter", jitter = 0.2)


##bee dry weight
par(mar = c(5,9,4,2))
boxplot(rbtw ~ species, data = t2, col = rainbow(16, alpha = 0.3),
        xlab = "weight (mg)", ylab = "",
        horizontal = TRUE, las = 1, range = 0, main = "Dry Weight")
stripchart(rbtw ~ species, data = t2, col = rainbow(16),
           pch = 19, cex = 2, las = 1, add = TRUE,
           method = "jitter", jitter = 0.2)
