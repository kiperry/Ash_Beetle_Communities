###########################################################
#
# Analysis of Ash Common Garden Beetle Data
#
# Focus on native herbivorous beetles
#
# NMDS analyses of bark and spray beetle data
# to compare communities on ash species and hybrid ash
#
# Rarefaction analyses of herbivore species by tree species
#
# Partition beta-diversity for bark and spray beetle data
#
# KI Perry; 11 June 2020; updated 15 June 2020
#
##########################################################

# install and load packages needed for NMDS and beta-diversity analyses

if (!suppressWarnings(require(vegan))) install.packages("vegan")
citation("vegan")

#install.packages("devtools")
#library(devtools)
#install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)
citation("pairwiseAdonis")

if (!suppressWarnings(require(betapart))) install.packages("betapart")
citation("betapart")

if (!suppressWarnings(require(BiodiversityR))) install.packages("BiodiversityR")
citation("BiodiversityR")

if (!suppressWarnings(require(ggplot2))) install.packages("ggplot2")
citation("ggplot2")

if (!suppressWarnings(require(reshape2))) install.packages("reshape2")
citation("reshape2")

if (!suppressWarnings(require(viridis))) install.packages("viridis")
citation("viridis")

# import data 

dat <- read.csv("./Data_Bark_Spray_NHerbivores.csv", row.names = 1)

# check structure of data

str(dat)
head(dat)
tail(dat)
dim(dat)

# change Block from an integer to a factor

dat <- within(dat, {
  block <- factor(block)
})
str(dat)

# set up color and point vectors for figures
colvec1 <- c("darkorange1", "darkgoldenrod1", "black", "chartreuse4")
colvec2 <- c("#FDE725FF", "#2D708EFF", "#481567FF", "#73D055FF")
pchvec1 <- c(16, 17, 15, 18)
pchvec2 <- c(19, 15, 4, 9)
ltyvec <- c(1, 2, 3, 4)

#########################################################

## Comparison of BARK wrap beetle data

# pull out bark data

bark <- dat[which(dat$method == "Bark"),]
str(bark)

#write.csv(bark, file = "BarkHBeetleData.csv")

# remove columns with zero values
# removes species not collected on bark wraps

bark <- bark[, colSums(bark != 0) > 0]

## Compare beetle species richness across tree species

# individual-based rarefaction by tree species, jackknife estimates by tree species
# first-order jackknife estimates are based on the number of singletons
# second-order jackknife estmiates are based on the number of singletons and doubletons

levels(bark$species)
rare.bark.H <- bark[which(bark$species == "H"),]
rare.bark.M <- bark[which(bark$species == "M"),]
rare.bark.N <- bark[which(bark$species == "N"),]
rare.bark.P <- bark[which(bark$species == "P"),]

sp.bark.H <- specaccum(rare.bark.H[7:64], method = "rarefaction", permutations = 100, gamma = "jack2")
sp.bark.M <- specaccum(rare.bark.M[7:64], method = "rarefaction", permutations = 100, gamma = "jack2")
sp.bark.N <- specaccum(rare.bark.N[7:64], method = "rarefaction", permutations = 100, gamma = "jack2")
sp.bark.P <- specaccum(rare.bark.P[7:64], method = "rarefaction", permutations = 100, gamma = "jack2")

plot(sp.bark.H, pch = 19, col = "#FDE725FF", xvar = c("individuals"), lty = 1, lwd = 3,
          ylab = "Species Richness", xlab = "Number of Individuals", xlim = c(0, 80), ylim = c(0, 35))
plot(sp.bark.M, add = TRUE, pch = 15, xvar = c("individuals"), lty = 2, lwd = 3, col = "#2D708EFF")
plot(sp.bark.N, add = TRUE, pch = 4, xvar = c("individuals"), lty = 3, lwd = 3, col = "#481567FF")
plot(sp.bark.P, add = TRUE, pch = 9, xvar = c("individuals"), lty = 4, lwd = 3, col = "#73D055FF")
legend("topleft", legend = c("BxM Hybrid","F. mandshurica","F. nigra", "F. pennsylvanica"), 
       lty = ltyvec, bty = "n", lwd = 3, col = colvec2)

levels(bark$species)

#calculates species richness for each site
specnumber(bark[7:64])

#calculates species richness by treatment
specnumber(bark[7:64], groups = bark$species)

bark.sp.j1 <- diversitycomp(bark[7:64], y = bark, factor1 = "species", index = "jack1")
bark.sp.j1
bark.sp.j2 <- diversitycomp(bark[7:64], y = bark, factor1 = "species", index = "jack2")
bark.sp.j2

bark.j1 <- diversityresult(bark[7:64], y=NULL, index = "jack1")
bark.j1
bark.j2 <- diversityresult(bark[7:64], y=NULL, index = "jack2")
bark.j2


#########################################################

## Partition beta-diversity for herbivores found on bark wraps

# pool abundance data across months for each block
# then combine all blocks into one dataset for analysis
# include appropriate site codes as row names

bark.B1 <- bark[which(bark$block == "1"),]
bark.B1.1 <- bark.B1[,-1:-6]
B_Block1 <- rowsum.data.frame(bark.B1.1, bark.B1$species)
B_Block1$code <- c("B1H", "B1M", "B1N", "B1P")
rownames(B_Block1) <- B_Block1[,59]
B_Block1 <- B_Block1[,-59]
rm(bark.B1, bark.B1.1)

bark.B2 <- bark[which(bark$block == "2"),]
bark.B2.2 <- bark.B2[,-1:-6]
B_Block2 <- rowsum.data.frame(bark.B2.2, bark.B2$species)
B_Block2$code <- c("B2H", "B2M", "B2N", "B2P")
rownames(B_Block2) <- B_Block2[,59]
B_Block2 <- B_Block2[,-59]
rm(bark.B2, bark.B2.2)

bark.B3 <- bark[which(bark$block == "3"),]
bark.B3.3 <- bark.B3[,-1:-6]
B_Block3 <- rowsum.data.frame(bark.B3.3, bark.B3$species)
B_Block3$code <- c("B3H", "B3M", "B3N", "B3P")
rownames(B_Block3) <- B_Block3[,59]
B_Block3 <- B_Block3[,-59]
rm(bark.B3, bark.B3.3)

bark.B4 <- bark[which(bark$block == "4"),]
bark.B4.4 <- bark.B4[,-1:-6]
B_Block4 <- rowsum.data.frame(bark.B4.4, bark.B4$species)
B_Block4$code <- c("B4H", "B4M", "B4N", "B4P")
rownames(B_Block4) <- B_Block4[,59]
B_Block4 <- B_Block4[,-59]
rm(bark.B4, bark.B4.4)

bark.B5 <- bark[which(bark$block == "5"),]
bark.B5.5 <- bark.B5[,-1:-6]
B_Block5 <- rowsum.data.frame(bark.B5.5, bark.B5$species)
B_Block5$code <- c("B5H", "B5M", "B5N", "B5P")
rownames(B_Block5) <- B_Block5[,59]
B_Block5 <- B_Block5[,-59]
rm(bark.B5, bark.B5.5)

bark2 <-rbind(B_Block1, B_Block2, B_Block3, B_Block4, B_Block5)
rm(B_Block1, B_Block2, B_Block3, B_Block4, B_Block5)

# convert matrix to presence/absence
bark2[bark2>0]<- 1

# add in grouping variables
bark2$species <- c("H", "M", "N", "P",
                   "H", "M", "N", "P",
                   "H", "M", "N", "P",
                   "H", "M", "N", "P",
                   "H", "M", "N", "P")

bark2$block <- c("1", "1", "1", "1",
                 "2", "2", "2", "2",
                 "3", "3", "3", "3",
                 "4", "4", "4", "4",
                 "5", "5", "5", "5")

bark2 <- bark2[,c(59,60,1:58)]

str(bark2)
bark2 <- within(bark2, {
  species <- factor(species)
  block <- factor(block)
})

# create beta part object for analyses
bark.core <- betapart.core(bark2[3:60])

# returns three dissimilarity matrices containing 
# pairwise between-site values of each beta-diversity component
bark2.dist <- beta.pair(bark.core, index.family = "sorensen")
str(bark2.dist)

# run NMDS models for each beta-diversity component

## Beta.sor - Total beta-diversity
nmds.bark.sor <- metaMDS(bark2.dist$beta.sor, trymax = 500, autotransform = TRUE)
nmds.bark.sor
stressplot(nmds.bark.sor)
goodness(nmds.bark.sor)
nmds.bark.sor$stress #stress is quality of fit
plot(nmds.bark.sor)

colors()
levels(bark2$species)
# H = black ash x Manchurian ash hybrid
# M = Manchurian ash (Fraxinus mandshurica)
# N = black ash 'Fall Gold' (Fraxinus nigra)
# P = green ash 'Patmore' (Fraxinus pennsylvanica)


ordiplot(nmds.bark.sor, type="n", xlim = c(-1, 0.8), ylim = c(-0.6, 0.6))
points(nmds.bark.sor, display = "sites", pch = pchvec1[bark2$species], cex = 1.5, 
       col = colvec2[bark2$species])
ordiellipse(nmds.bark.sor, groups = bark2$species, display = "sites", draw = "lines", 
            col = colvec2[bark2$species], lwd = 3, conf = 0.90)
ordihull(nmds.bark.sor, groups = bark2$species, display = "sites", draw = c("polygon"), col = NULL,
         border = colvec2[bark2$species], lwd = 2.5)
legend("topleft", legend = c("BxM hybrid","F. mandshurica","F. nigra", "F. pennsylvanica"), 
       pch = pchvec1[bark2$species], cex = 1.2, bty = "n", col = colvec2[bark2$species])


## Test for differences in herbivore beetle composition across ash tree species
## for tota beta-diversity component

# PERMANOVA tests whether the group centroid of beetle communities on tree species
# differs in multivariate space (e.g. different community composition)
adonis2(bark2.dist$beta.sor ~ bark2$species, permutations = 999)
adonis2(bark2.dist$beta.sor ~ bark2$block, permutations = 999)
pairwise.adonis(bark2.dist$beta.sor, bark2$species)

# BETADISPER tests whether the dispersion of a treatment from its spatial median is different
# between groups (i.e. species redundancy across space)
# analysis of multivariate homogeneity of group dispersions (variances)
# multivariate analogue of Levene's test for homogenetiy of variances
bark.beta.sor <- betadisper(bark2.dist$beta.sor, bark2$species, type = c("median"))
bark.beta.sor
anova(bark.beta.sor)
plot(bark.beta.sor)
boxplot(bark.beta.sor, ylab = "Distance to median")
TukeyHSD(bark.beta.sor, which = "group", conf.level = 0.95)


## Beta.sim - Turnover
nmds.bark.sim <- metaMDS(bark2.dist$beta.sim, trymax = 500, autotransform = TRUE)
nmds.bark.sim
stressplot(nmds.bark.sim)
goodness(nmds.bark.sim)
nmds.bark.sim$stress
plot(nmds.bark.sim) 

ordiplot(nmds.bark.sim, type="n", xlim = c(-1, 0.8), ylim = c(-0.6, 0.6))
points(nmds.bark.sim, display = "sites", pch = pchvec1[bark2$species], cex = 1.5, 
       col = colvec2[bark2$species])
ordiellipse(nmds.bark.sim, groups = bark2$species, display = "sites", draw = "lines", 
            col = colvec2[bark2$species], lwd = 3, conf = 0.90)
ordihull(nmds.bark.sim, groups = bark2$species, display = "sites", draw = c("polygon"), col = NULL,
         border = colvec2[bark2$species], lwd = 2.5)
legend("topright", legend = c("BxM Hybrid","F. mandshurica","F. nigra", "F. pennsylvanica"), 
       pch = pchvec1[bark2$species], cex = 1.2, bty = "n", col = colvec2[bark2$species])


## Test for differences in herbivore beetle composition across ash tree species
## for turnover beta-diversity component

# PERMANOVA tests whether the group centroid of beetle communities on tree species
# differs in multivariate space (e.g. different community composition)
adonis2(bark2.dist$beta.sim ~ bark2$species, permutations = 999)
adonis2(bark2.dist$beta.sim ~ bark2$block, permutations = 999)
pairwise.adonis(bark2.dist$beta.sim, bark2$species)

# BETADISPER tests whether the dispersion of a treatment from its spatial median is different
# between groups
# analysis of multivariate homogeneity of group dispersions (variances)
# multivariate analogue of Levene's test for homogenetiy of variances
bark.beta.sim <- betadisper(bark2.dist$beta.sim, bark2$species, type = c("median"))
bark.beta.sim
anova(bark.beta.sim)
plot(bark.beta.sim)
boxplot(bark.beta.sim, ylab = "Distance to median")
TukeyHSD(bark.beta.sim, which = "group", conf.level = 0.95)


## Beta.sne - Nestedness
nmds.bark.sne <- metaMDS(bark2.dist$beta.sne, trymax = 500, autotransform = TRUE)
nmds.bark.sne
stressplot(nmds.bark.sne)
goodness(nmds.bark.sne)
nmds.bark.sne$stress
plot(nmds.bark.sne)

ordiplot(nmds.bark.sne, type="n", xlim = c(-0.4, 0.4), ylim = c(-0.6, 0.4))
points(nmds.bark.sne, display = "sites", pch = pchvec1[bark2$species], cex = 1.5, 
       col = colvec2[bark2$species])
ordiellipse(nmds.bark.sne, groups = bark2$species, display = "sites", draw = "lines", 
            col = colvec2[bark2$species], lwd = 3, conf = 0.90)
ordihull(nmds.bark.sne, groups = bark2$species, display = "sites", draw = c("polygon"), col = NULL,
         border = colvec2[bark2$species], lwd = 2.5)
legend("topright", legend = c("BxM Hybrid","F. mandshurica","F. nigra", "F. pennsylvanica"), 
       pch = pchvec[bark2$species], cex = 1.2, bty = "n", col = colvec[bark2$species])


## Test for differences in herbivore beetle composition across ash tree species
## for nestedness beta-diversity component

# PERMANOVA tests whether the group centroid of beetle communities on tree species
# differs in multivariate space (e.g. different community composition)
adonis2(bark2.dist$beta.sne ~ bark2$species, permutations = 999)
adonis2(bark2.dist$beta.sne ~ bark2$block, permutations = 999)
pairwise.adonis(bark2.dist$beta.sne, bark2$species)

# BETADISPER tests whether the dispersion of a treatment from its spatial median is different
# between groups (i.e. species redundancy across space)
# analysis of multivariate homogeneity of group dispersions (variances)
# multivariate analogue of Levene's test for homogenetiy of variances
bark.beta.sne <- betadisper(bark2.dist$beta.sne, bark2$species, type = c("median"))
bark.beta.sne
anova(bark.beta.sne)
plot(bark.beta.sne)
boxplot(bark.beta.sne, ylab = "Distance to median")
TukeyHSD(bark.beta.sne, which = "group", conf.level = 0.95)


#########################################################
## Descriptive total dissimilarity across all trees by species
# returns the total dissimilarity across all n sites
# along with turnover and nestedness components of that dissimilarity

# pull out reps of each trees species
bark.N <- bark2[which(bark2$species == "N"),]
bark.H <- bark2[which(bark2$species == "H"),]
bark.M <- bark2[which(bark2$species == "M"),]
bark.P <- bark2[which(bark2$species == "P"),]

bark.core.N <- betapart.core(bark.N[3:60])
bark.core.H <- betapart.core(bark.H[3:60])
bark.core.M <- betapart.core(bark.M[3:60])
bark.core.P <- betapart.core(bark.P[3:60])

beta.multi(bark.core.H, index.family = "sorensen")
beta.multi(bark.core.M, index.family = "sorensen")
beta.multi(bark.core.N, index.family = "sorensen")
beta.multi(bark.core.P, index.family = "sorensen")

# create stacked bar plots for each species

theme_set(theme_bw())

# create a data frame that contains the beta.sim and beta.sne
# for each ash species from the beta.multi outputs above
# beta.sim = species turnover; beta.sne = species nestedness
# beta.sor = total betadiversity (i.e., beta.sim + beta.sne = beta.sor)

# bark beta-diversity
b.b <- data.frame(
  Turnover = c(0.78, 0.68, 0.73, 0.75),
  Nestedness = c(0.05, 0.12, 0.06, 0.04)
)

b.b$species <- c("NxM hybrid","F. mandshurica","F. nigra", "F. pennsylvanica")

b.b <- melt(b.b, id.vars = "species")

png("NBetaDiversity_Bark_Ash.png", width = 1000, height = 800)
ggplot(b.b, aes(fill = variable, x = species, y = value)) + 
  geom_bar(position = "stack", stat = "identity") +
  ylab("Beta-diversity") + xlab("") +
  coord_flip() + theme(text = element_text(size = 20)) +
  theme(legend.title = element_blank()) + ylim(0, max(1)) + 
  scale_fill_manual(values=c("gray64", "gray14"))
dev.off()


#########################################################

# Comparison of SPRAY beetle data

# pull out spray data

spray <- dat[which(dat$method == "Spray"),]
str(spray)

# remove columns with zero values
# removes species not collected via spraying

spray <- spray[, colSums(spray != 0) > 0]

## Compare beetle species richness across tree species

# individual-based rarefaction by tree species, jackknife estimates by tree species
# first-order jackknife estimates are based on the number of singletons
# second-order jackknife estmiates are based on the number of singletons and doubletons

levels(spray$species)
rare.spray.H <- spray[which(spray$species == "H"),]
rare.spray.M <- spray[which(spray$species == "M"),]
rare.spray.N <- spray[which(spray$species == "N"),]
rare.spray.P <- spray[which(spray$species == "P"),]

sp.spray.H <- specaccum(rare.spray.H[7:45], method = "rarefaction", permutations = 100, gamma = "jack2")
sp.spray.M <- specaccum(rare.spray.M[7:45], method = "rarefaction", permutations = 100, gamma = "jack2")
sp.spray.N <- specaccum(rare.spray.N[7:45], method = "rarefaction", permutations = 100, gamma = "jack2")
sp.spray.P <- specaccum(rare.spray.P[7:45], method = "rarefaction", permutations = 100, gamma = "jack2")

plot(sp.spray.H, col = "#FDE725FF", xvar = c("individuals"), lty = 1, lwd = 3,
     ylab = "Species Richness", xlab = "Number of Individuals", xlim = c(0, 80), ylim = c(0, 35))
plot(sp.spray.M, add = TRUE, xvar = c("individuals"), lty = 2, lwd = 3, col = "#2D708EFF")
plot(sp.spray.N, add = TRUE, xvar = c("individuals"), lty = 3, lwd = 3, col = "#481567FF")
plot(sp.spray.P, add = TRUE, xvar = c("individuals"), lty = 4, lwd = 3, col = "#73D055FF")
legend("topleft", legend = c("BxM Hybrid","F. mandshurica","F. nigra", "F. pennsylvanica"), 
       lty = ltyvec, bty = "n", lwd = 3, col = colvec2)

#calculates species richness for each site
specnumber(spray[7:45])

#calculates species richness by treatment
specnumber(spray[7:45], groups = spray$species)

spray.sp.j1 <- diversitycomp(spray[7:45], y = spray, factor1 = "species", index = "jack1")
spray.sp.j1
spray.sp.j2 <- diversitycomp(spray[7:45], y = spray, factor1 = "species", index = "jack2")
spray.sp.j2

spray.j1 <- diversityresult(spray[7:45], y=NULL, index = "jack1")
spray.j1
spray.j2 <- diversityresult(spray[7:45], y=NULL, index = "jack2")
spray.j2


#########################################################

## Partition beta-diversity for herbivores collected via spraying

# pool abundance data across months for each block
# then combine all blocks into one dataset for analysis
# include appropriate site codes as row names

spray.B1 <- spray[which(spray$block == "1"),]
spray.B1.1 <- spray.B1[,-1:-6]
S_Block1 <- rowsum.data.frame(spray.B1.1, spray.B1$species)
S_Block1$code <- c("B1H", "B1M", "B1N", "B1P")
rownames(S_Block1) <- S_Block1[,40]
S_Block1 <- S_Block1[,-40]
rm(spray.B1, spray.B1.1)

spray.B2 <- spray[which(spray$block == "2"),]
spray.B2.2 <- spray.B2[,-1:-6]
S_Block2 <- rowsum.data.frame(spray.B2.2, spray.B2$species)
S_Block2$code <- c("B2H", "B2M", "B2N", "B2P")
rownames(S_Block2) <- S_Block2[,40]
S_Block2 <- S_Block2[,-40]
rm(spray.B2, spray.B2.2)

spray.B3 <- spray[which(spray$block == "3"),]
spray.B3.3 <- spray.B3[,-1:-6]
S_Block3 <- rowsum.data.frame(spray.B3.3, spray.B3$species)
S_Block3$code <- c("B3H", "B3M", "B3N", "B3P")
rownames(S_Block3) <- S_Block3[,40]
S_Block3 <- S_Block3[,-40]
rm(spray.B3, spray.B3.3)

spray.B4 <- spray[which(spray$block == "4"),]
spray.B4.4 <- spray.B4[,-1:-6]
S_Block4 <- rowsum.data.frame(spray.B4.4, spray.B4$species)
S_Block4$code <- c("B4H", "B4M", "B4N", "B4P")
rownames(S_Block4) <- S_Block4[,40]
S_Block4 <- S_Block4[,-40]
rm(spray.B4, spray.B4.4)

spray.B5 <- spray[which(spray$block == "5"),]
spray.B5.5 <- spray.B5[,-1:-6]
S_Block5 <- rowsum.data.frame(spray.B5.5, spray.B5$species)
S_Block5$code <- c("B5H", "B5M", "B5N", "B5P")
rownames(S_Block5) <- S_Block5[,40]
S_Block5 <- S_Block5[,-40]
rm(spray.B5, spray.B5.5)

spray2 <-rbind(S_Block1, S_Block2, S_Block3, S_Block4, S_Block5)
rm(S_Block1, S_Block2, S_Block3, S_Block4, S_Block5)

# convert matrix to presence/absence
spray2[spray2>0]<- 1

# add in grouping variables
spray2$species <- c("H", "M", "N", "P",
                   "H", "M", "N", "P",
                   "H", "M", "N", "P",
                   "H", "M", "N", "P",
                   "H", "M", "N", "P")

spray2$block <- c("1", "1", "1", "1",
                 "2", "2", "2", "2",
                 "3", "3", "3", "3",
                 "4", "4", "4", "4",
                 "5", "5", "5", "5")

spray2 <- spray2[,c(40,41,1:39)]

str(spray2)
spray2 <- within(spray2, {
  species <- factor(species)
  block <- factor(block)
})

# create beta part object for analyses
spray.core <- betapart.core(spray2[3:41])

# returns three dissimilarity matrices containing 
# pairwise between-site values of each beta-diversity component
spray2.dist <- beta.pair(spray.core, index.family = "sorensen")
str(spray2.dist)

# run NMDS models for each beta-diversity component

## Beta.sor - Total beta-diversity
nmds.spray.sor <- metaMDS(spray2.dist$beta.sor, trymax = 500, autotransform = TRUE)
nmds.spray.sor
stressplot(nmds.spray.sor)
goodness(nmds.spray.sor)
nmds.spray.sor$stress #stress is quality of fit
plot(nmds.spray.sor)

colors()
levels(spray2$species)
# H = black ash x Manchurian ash hybrid
# M = Manchurian ash (Fraxinus mandshurica)
# N = black ash 'Fall Gold' (Fraxinus nigra)
# P = green ash 'Patmore' (Fraxinus pennsylvanica)


ordiplot(nmds.spray.sor, type="n", xlim = c(-1.0, 1.0), ylim = c(-0.8, 0.8))
points(nmds.spray.sor, display = "sites", pch = pchvec1[spray2$species], cex = 1.5, 
       col = colvec2[spray2$species])
ordiellipse(nmds.spray.sor, groups = spray2$species, display = "sites", draw = "lines", 
            col = colvec2[spray2$species], lwd = 2.5, conf = 0.90)
ordihull(nmds.spray.sor, groups = spray2$species, display = "sites", draw = c("polygon"), col = NULL,
         border = colvec2[spray2$species], lwd = 2.5)
legend("topleft", legend = c("BxM Hybrid","F. mandshurica","F. nigra", "F. pennsylvanica"), 
       pch = pchvec1[spray2$species], cex = 1.2, bty = "n", col = colvec2[spray2$species])


## Test for differences in herbivore beetle composition across ash tree species
## for total beta-diversity component

# PERMANOVA tests whether the group centroid of beetle communities on tree species
# differs in multivariate space (e.g. different community composition)
adonis2(spray2.dist$beta.sor ~ spray2$species, permutations = 999)
adonis2(spray2.dist$beta.sor ~ spray2$block, permutations = 999)
pairwise.adonis(spray2.dist$beta.sor, spray2$species)

# BETADISPER tests whether the dispersion of a treatment from its spatial median is different
# between groups (i.e. species redundancy across space)
# analysis of multivariate homogeneity of group dispersions (variances)
# multivariate analogue of Levene's test for homogenetiy of variances
spray.beta.sor <- betadisper(spray2.dist$beta.sor, spray2$species, type = c("median"))
spray.beta.sor
anova(spray.beta.sor)
plot(spray.beta.sor)
boxplot(spray.beta.sor, ylab = "Distance to median")
TukeyHSD(spray.beta.sor, which = "group", conf.level = 0.95)


## Beta.sim - Turnover
nmds.spray.sim <- metaMDS(spray2.dist$beta.sim, trymax = 500, autotransform = TRUE)
nmds.spray.sim
stressplot(nmds.spray.sim)
goodness(nmds.spray.sim)
nmds.spray.sim$stress
plot(nmds.spray.sim) 

ordiplot(nmds.spray.sim, type="n", xlim = c(-1, 1), ylim = c(-0.8, 0.8))
points(nmds.spray.sim, display = "sites", pch = pchvec1[spray2$species], cex = 1.5, 
       col = colvec2[spray2$species])
ordiellipse(nmds.spray.sim, groups = spray2$species, display = "sites", draw = "lines", 
            col = colvec2[spray2$species], lwd = 2.5, conf = 0.90)
ordihull(nmds.spray.sim, groups = spray2$species, display = "sites", draw = c("polygon"), col = NULL,
         border = colvec2[spray2$species], lwd = 2.5)
legend("topright", legend = c("BxM Hybrid","F. mandshurica","F. nigra", "F. pennsylvanica"), 
       pch = pchvec[spray2$species], cex = 1.2, bty = "n", col = colvec[spray2$species])


## Test for differences in herbivore beetle composition across ash tree species
## for turnover beta-diversity component

# PERMANOVA tests whether the group centroid of beetle communities on tree species
# differs in multivariate space (e.g. different community composition)
adonis2(spray2.dist$beta.sim ~ spray2$species, permutations = 999)
adonis2(spray2.dist$beta.sim ~ spray2$block, permutations = 999)
pairwise.adonis(spray2.dist$beta.sim, spray2$species)

# BETADISPER tests whether the dispersion of a treatment from its spatial median is different
# between groups
# analysis of multivariate homogeneity of group dispersions (variances)
# multivariate analogue of Levene's test for homogenetiy of variances
spray.beta.sim <- betadisper(spray2.dist$beta.sim, spray2$species, type = c("median"))
spray.beta.sim
anova(spray.beta.sim)
plot(spray.beta.sim)
boxplot(spray.beta.sim, ylab = "Distance to median")
TukeyHSD(spray.beta.sim, which = "group", conf.level = 0.95)


## Beta.sne - Nestedness
nmds.spray.sne <- metaMDS(spray2.dist$beta.sne, trymax = 500, autotransform = TRUE)
nmds.spray.sne
stressplot(nmds.spray.sne)
goodness(nmds.spray.sne)
nmds.spray.sne$stress
plot(nmds.spray.sne)

ordiplot(nmds.spray.sne, type="n", xlim = c(-0.3, 0.3), ylim = c(-0.3, 0.3))
points(nmds.spray.sne, display = "sites", pch = pchvec1[spray2$species], cex = 1.5, 
       col = colvec2[spray2$species])
ordiellipse(nmds.spray.sne, groups = spray2$species, display = "sites", draw = "lines", 
            col = colvec2[spray2$species], lwd = 2.5, conf = 0.90)
ordihull(nmds.spray.sne, groups = spray2$species, display = "sites", draw = c("polygon"), col = NULL,
         border = colvec2[spray2$species], lwd = 2.5)
legend("topright", legend = c("BxM Hybrid","F. mandshurica","F. nigra", "F. pennsylvanica"), 
       pch = pchvec[spray2$species], cex = 1.2, bty = "n", col = colvec[spray2$species])


## Test for differences in herbivore beetle composition across ash tree species
## for nestedness beta-diversity component

# PERMANOVA tests whether the group centroid of beetle communities on tree species
# differs in multivariate space (e.g. different community composition)
adonis2(spray2.dist$beta.sne ~ spray2$species, permutations = 999)
adonis2(spray2.dist$beta.sne ~ spray2$block, permutations = 999)
pairwise.adonis(spray2.dist$beta.sne, spray2$species)

# BETADISPER tests whether the dispersion of a treatment from its spatial median is different
# between groups
# analysis of multivariate homogeneity of group dispersions (variances)
# multivariate analogue of Levene's test for homogenetiy of variances
spray.beta.sne <- betadisper(spray2.dist$beta.sne, spray2$species, type = c("median"))
spray.beta.sne
anova(spray.beta.sne)
plot(spray.beta.sne)
boxplot(spray.beta.sne, ylab = "Distance to median")
TukeyHSD(spray.beta.sne, which = "group", conf.level = 0.95)


#########################################################

## Descriptive total dissimilarity across all trees by species
# returns the total dissimilarity across all n sites
# along with turnover and nestedness components of that dissimilarity

# pull out reps of each trees species
spray.H <- spray2[which(spray2$species == "H"),]
spray.M <- spray2[which(spray2$species == "M"),]
spray.N <- spray2[which(spray2$species == "N"),]
spray.P <- spray2[which(spray2$species == "P"),]

spray.core.H <- betapart.core(spray.H[3:41])
spray.core.M <- betapart.core(spray.M[3:41])
spray.core.N <- betapart.core(spray.N[3:41])
spray.core.P <- betapart.core(spray.P[3:41])

beta.multi(spray.core.H, index.family = "sorensen")
beta.multi(spray.core.M, index.family = "sorensen")
beta.multi(spray.core.N, index.family = "sorensen")
beta.multi(spray.core.P, index.family = "sorensen")

# create a data frame that contains the beta.sim and beta.sne
# for each ash species from the beta.multi outputs above
# beta.sim = species turnover; beta.sne = species nestedness
# beta.sor = total betadiversity (i.e., beta.sim + beta.sne = beta.sor)

# spray beta-diversity
s.b <- data.frame(
  Turnover = c(0.65, 0.64, 0.61, 0.72),
  Nestedness = c(0.08, 0.08, 0.07, 0.06)
)

s.b$species <- c("NxM Hybrid","F. mandshurica","F. nigra", "F. pennsylvanica")

s.b <- melt(s.b, id.vars = "species")

png("NBetaDiversity_Spray_Ash.png", width = 1000, height = 800)
ggplot(s.b, aes(fill = variable, x = species, y = value)) + 
  geom_bar(position = "stack", stat = "identity") +
  ylab("Beta-diversity") + xlab("") +
  coord_flip() + theme(text = element_text(size = 20)) +
  theme(legend.title = element_blank()) + ylim(0, max(1)) + 
  scale_fill_manual(values=c("gray64", "gray14"))
dev.off()


###############################################################################################

## Make rarefaction figure
png("NSpecies_Rarefaction_Ash.png", width = 2000, height = 800)
layout(matrix(c(1,2), nrow = 1, ncol = 2, byrow = TRUE))
layout.show(2)

par(mar=c(6,7,3,1))
plot(sp.spray.H, col = "#FDE725FF", xvar = c("individuals"), lty = 1, lwd = 3, cex.lab = 2.5, cex.axis = 1.8,
     ylab = "Species Richness", xlab = "Number of Individuals", xlim = c(0, 80), ylim = c(0, 35))
plot(sp.spray.M, add = TRUE, xvar = c("individuals"), lty = 2, lwd = 3, col = "#2D708EFF")
plot(sp.spray.N, add = TRUE, xvar = c("individuals"), lty = 3, lwd = 3, col = "#481567FF")
plot(sp.spray.P, add = TRUE, xvar = c("individuals"), lty = 4, lwd = 3, col = "#73D055FF")
legend("topleft", legend = c("BxM hybrid","F. mandshurica","F. nigra", "F. pennsylvanica"), 
       lty = ltyvec, bty = "n", lwd = 5, col = colvec2, cex = 2.5)
text(76,34, "A", pos = 4, font = 2, cex = 3)

par(mar=c(6,7,3,1))
plot(sp.bark.H, col = "#FDE725FF", xvar = c("individuals"), lty = 1, lwd = 3, cex.lab = 2.5, cex.axis = 1.8,
     ylab = "Species Richness", xlab = "Number of Individuals", xlim = c(0, 80), ylim = c(0, 35))
plot(sp.bark.M, add = TRUE, xvar = c("individuals"), lty = 2, lwd = 3, col = "#2D708EFF")
plot(sp.bark.N, add = TRUE, xvar = c("individuals"), lty = 3, lwd = 3, col = "#481567FF")
plot(sp.bark.P, add = TRUE, xvar = c("individuals"), lty = 4, lwd = 3, col = "#73D055FF")
text(76,34, "B", pos = 4, font = 2, cex = 3)
dev.off()


################################################################################################

## Make spray NMDS figure

png("NMDS_Herbivores_Spray_Ash.png", width = 800, height = 2000, pointsize = 30)
layout(matrix(c(1,2,3), nrow = 3, ncol = 1, byrow = TRUE))
layout.show(3)

par(mar=c(5,4,3,2))
ordiplot(nmds.spray.sor, type="n", xlim = c(-1.5, 1.0), ylim = c(-1, 1))
points(nmds.spray.sor, display = "sites", pch = pchvec1[spray2$species], cex = 1.5, 
       col = colvec2[spray2$species])
ordiellipse(nmds.spray.sor, groups = spray2$species, display = "sites", draw = "lines", 
            col = colvec2[spray2$species], lwd = 3, conf = 0.90)
legend("topleft", legend = c("BxM hybrid","F. mandshurica","F. nigra", "F. pennsylvanica"), 
       pch = pchvec1[spray2$species], cex = 1.5, bty = "n", col = colvec2[spray2$species])
text(1, 0.85, "A", pos = 4, font = 2, cex = 1.8)

ordiplot(nmds.spray.sim, type="n", xlim = c(-1, 1), ylim = c(-0.8, 0.8))
points(nmds.spray.sim, display = "sites", pch = pchvec1[spray2$species], cex = 1.5, 
       col = colvec2[spray2$species])
ordiellipse(nmds.spray.sim, groups = spray2$species, display = "sites", draw = "lines", 
            col = colvec2[spray2$species], lwd = 3, conf = 0.90)
text(1, 0.70, "B", pos = 4, font = 2, cex = 1.8)

ordiplot(nmds.spray.sne, type="n", xlim = c(-0.3, 0.3), ylim = c(-0.3, 0.3))
points(nmds.spray.sne, display = "sites", pch = pchvec1[spray2$species], cex = 1.5, 
       col = colvec2[spray2$species])
ordiellipse(nmds.spray.sne, groups = spray2$species, display = "sites", draw = "lines", 
            col = colvec2[spray2$species], lwd = 3, conf = 0.90)
text(0.375, 0.265, "C", pos = 4, font = 2, cex = 1.8)

dev.off()


############################################################################

# Make bark NMDS figure

png("NMDS_Herbivores_Bark_Ash.png", width = 800, height = 2000, pointsize = 30)
layout(matrix(c(1,2,3), nrow = 3, ncol = 1, byrow = TRUE))
layout.show(3)

par(mar=c(5,4,3,2))

ordiplot(nmds.bark.sor, type="n", xlim = c(-1.5, 0.8), ylim = c(-0.6, 0.6))
points(nmds.bark.sor, display = "sites", pch = pchvec1[bark2$species], cex = 1.5, 
       col = colvec2[bark2$species])
ordiellipse(nmds.bark.sor, groups = bark2$species, display = "sites", draw = "lines", 
            col = colvec2[bark2$species], lwd = 3, conf = 0.90)
legend("topleft", legend = c("BxM hybrid","F. mandshurica","F. nigra", "F. pennsylvanica"), 
       pch = pchvec1[bark2$species], cex = 1.5, bty = "n", col = colvec2[bark2$species])
text(0.66, 0.75, "A", pos = 4, font = 2, cex = 1.8)

ordiplot(nmds.bark.sim, type="n", xlim = c(-1, 0.8), ylim = c(-0.6, 0.6))
points(nmds.bark.sim, display = "sites", pch = pchvec1[bark2$species], cex = 1.5, 
       col = colvec2[bark2$species])
ordiellipse(nmds.bark.sim, groups = bark2$species, display = "sites", draw = "lines", 
            col = colvec2[bark2$species], lwd = 3, conf = 0.90)
text(0.68, 0.57, "B", pos = 4, font = 2, cex = 1.8)

ordiplot(nmds.bark.sne, type="n", xlim = c(-0.4, 0.4), ylim = c(-0.6, 0.4))
points(nmds.bark.sne, display = "sites", pch = pchvec1[bark2$species], cex = 1.5, 
       col = colvec2[bark2$species])
ordiellipse(nmds.bark.sne, groups = bark2$species, display = "sites", draw = "lines", 
            col = colvec2[bark2$species], lwd = 3, conf = 0.90)
text(0.60, 0.34, "C", pos = 4, font = 2, cex = 1.8)

dev.off()

#################################################################################

# Combined spray and bark NMDS figure

png("NMDS_NHerbivores_Canopy_Bark_Ash.png", width = 1600, height = 2000, pointsize = 30)

par(mfrow=c(3,2))
par(mar=c(5,4,3,2))

ordiplot(nmds.spray.sor, type="n", xlim = c(-1.5, 1.0), ylim = c(-1, 1))
points(nmds.spray.sor, display = "sites", pch = pchvec1[spray2$species], cex = 1.5, 
       col = colvec2[spray2$species])
ordiellipse(nmds.spray.sor, groups = spray2$species, display = "sites", draw = "lines", 
            col = colvec2[spray2$species], lwd = 3, conf = 0.90)
legend("topleft", legend = c("BxM hybrid","F. mandshurica","F. nigra", "F. pennsylvanica"), 
       pch = pchvec1[bark2$species], cex = 1.5, bty = "n", col = colvec2[bark2$species])
text(0.99, 0.89, "A", pos = 4, font = 2, cex = 1.8)

ordiplot(nmds.bark.sor, type="n", xlim = c(-1.5, 0.8), ylim = c(-0.6, 0.6))
points(nmds.bark.sor, display = "sites", pch = pchvec1[bark2$species], cex = 1.5, 
       col = colvec2[bark2$species])
ordiellipse(nmds.bark.sor, groups = bark2$species, display = "sites", draw = "lines", 
            col = colvec2[bark2$species], lwd = 3, conf = 0.90)
text(0.66, 0.75, "B", pos = 4, font = 2, cex = 1.8)

ordiplot(nmds.spray.sim, type="n", xlim = c(-1, 1), ylim = c(-0.8, 0.8))
points(nmds.spray.sim, display = "sites", pch = pchvec1[spray2$species], cex = 1.5, 
       col = colvec2[spray2$species])
ordiellipse(nmds.spray.sim, groups = spray2$species, display = "sites", draw = "lines", 
            col = colvec2[spray2$species], lwd = 3, conf = 0.90)
text(0.99, 0.70, "C", pos = 4, font = 2, cex = 1.8)

ordiplot(nmds.bark.sim, type="n", xlim = c(-1, 0.8), ylim = c(-0.6, 0.6))
points(nmds.bark.sim, display = "sites", pch = pchvec1[bark2$species], cex = 1.5, 
       col = colvec2[bark2$species])
ordiellipse(nmds.bark.sim, groups = bark2$species, display = "sites", draw = "lines", 
            col = colvec2[bark2$species], lwd = 3, conf = 0.90)
text(0.68, 0.57, "D", pos = 4, font = 2, cex = 1.8)

ordiplot(nmds.spray.sne, type="n", xlim = c(-0.3, 0.3), ylim = c(-0.3, 0.3))
points(nmds.spray.sne, display = "sites", pch = pchvec1[spray2$species], cex = 1.5, 
       col = colvec2[spray2$species])
ordiellipse(nmds.spray.sne, groups = spray2$species, display = "sites", draw = "lines", 
            col = colvec2[spray2$species], lwd = 3, conf = 0.90)
text(0.37, 0.265, "E", pos = 4, font = 2, cex = 1.8)

ordiplot(nmds.bark.sne, type="n", xlim = c(-0.4, 0.4), ylim = c(-0.6, 0.4))
points(nmds.bark.sne, display = "sites", pch = pchvec1[bark2$species], cex = 1.5, 
       col = colvec2[bark2$species])
ordiellipse(nmds.bark.sne, groups = bark2$species, display = "sites", draw = "lines", 
            col = colvec2[bark2$species], lwd = 3, conf = 0.90)
text(0.60, 0.34, "F", pos = 4, font = 2, cex = 1.8)

dev.off()