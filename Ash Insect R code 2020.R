#Packages
library(lme4)
library(lmerTest)
library(blmeco)
library(car)
library(multcomp)
library(pscl)

#Enter data
data <- read.csv("data_5.8.20.csv", header=TRUE)
data

#Set Period, Block, Tree and year as factors vs. integers
data$Period <- as.factor(data$Period)
data$Block <- as.factor(data$Block)
data$Tree <- as.factor(data$Tree)

#Chop up data
bark <- data[which(data$Method == "bk"),]
spray <- data[which(data$Method == "sp"),]

data$natabun <-(data$Tot_abun - data$Invasive)
data$natherbabun <-(data$Herb_abun - data$InvasHerb_abun)

#####Bark; response: Total abundance
hist(bark$Tot_abun)

bkta1 <- glmer(Tot_abun ~ Species + (1|Block) , family=poisson, data = bark)
summary(bkta1)
plot(bkta1)
dispersion_glmer(bkta1)
Anova(bkta1)
qqnorm(resid(bkta1))
qqline(resid(bkta1))


bkta2 <- glmer(Tot_abun ~ Species + Period +  (1|Block), family=poisson, data = bark)
summary(bkta2)
plot(bkta2)
dispersion_glmer(bkta2)
Anova(bkta2)

bkta3 <- glmer(Tot_abun ~ Species * Period +  (1|Block), family=poisson, data = bark)
summary(bkta3)
plot(bkta3)
dispersion_glmer(bkta3)
Anova(bkta3)

lsmeans(lm3, pairwise~Month*Species, adjust="tukey")

anova(bkta3,bkta2)

bkta4 <- glmer.nb(Tot_abun ~ Species + (1|Block) , data = bark)
summary(bkta4)
plot(bkta4)
dispersion_glmer(bkta4)
Anova(bkta4)
anova(bkta4)
qqnorm(resid(bkta4))
qqline(resid(bkta4))
lsmeans(bkta4, pairwise~Species, adjust="tukey")

anova(bkta4,bkta1)

#Best
bkta5 <- glmer.nb(Tot_abun ~ Species + Period +  (1|Block), data = bark)
summary(bkta5)
plot(bkta5)
dispersion_glmer(bkta5)
Anova(bkta5)

anova(bkta5,bkta4)

#Bad
bkta6 <- glmer.nb(Tot_abun ~ Species * Period +  (1|Block), data = bark)
summary(bkta6)
plot(bkta6)
dispersion_glmer(bkta6)
Anova(bkta6)

#contrasts
summary(glht(bkta5, linfct = mcp(Species = "Tukey")))
summary(glht(bkta5, linfct = mcp(Period = "Tukey")))

K1 <- glht(bkta5, linfct = mcp(Species = "Tukey"))$linfct
K2 <- glht(bkta5, linfct = mcp(Period = "Tukey"))$linfct

summary(glht(bkta5, linfct = rbind(K1, K2)))


######Bark; response: Herbivore Abundance
hist(bark$Herb_abun)

bkha1 <- glmer(Herb_abun ~ Species + (1|Block) , family=poisson, data = bark)
summary(bkha1)
plot(bkha1)
dispersion_glmer(bkha1)
Anova(bkha1)

bkha2 <- glmer(Herb_abun ~ Species + Period +  (1|Block), family=poisson, data = bark)
summary(bkha2)
plot(bkha2)
dispersion_glmer(bkha2)
Anova(bkha2)

anova(bkha2,bkha1)

bkha3 <- glmer(Herb_abun ~ Species * Period +  (1|Block), family=poisson, data = bark)
summary(bkha3)
plot(bkha3)
dispersion_glmer(bkha3)
Anova(bkha3)

anova(bkha3,bkha2)

bkha4 <- glmer.nb(Herb_abun ~ Species + (1|Block) , data = bark)
summary(bkha4)
plot(bkha4)
dispersion_glmer(bkha4)
Anova(bkha4)
lsmeans(bkha4, pairwise~Species, adjust="tukey")

anova(bkha4,bkha1)

#Best
bkha5 <- glmer.nb(Herb_abun ~ Species + Period +  (1|Block), data = bark)
summary(bkha5)
plot(bkha5)
dispersion_glmer(bkha5)
Anova(bkha5)

anova(bkha5,bkha4)

#Bad
bkha6 <- glmer.nb(Herb_abun ~ Species * Period +  (1|Block), data = bark)
summary(bkha6)
plot(bkha6)
dispersion_glmer(bkha6)
Anova(bkha6)

anova(bkha6, bkha5)

#contrasts
summary(glht(bkha5, linfct = mcp(Species = "Tukey")))
summary(glht(bkha5, linfct = mcp(Period = "Tukey")))

K3 <- glht(bkha5, linfct = mcp(Species = "Tukey"))$linfct
K4 <- glht(bkha5, linfct = mcp(Period = "Tukey"))$linfct

summary(glht(bkha5, linfct = rbind(K3, K4)))


#####Specialist Abundance; Bark; Only model 5
hist(spray$natabun)
data
bks1 <- glmer(NatHerb_Rich ~ Species + (1|Block) , family=poisson, data = spray)
summary(bks1)
plot(bks1)
dispersion_glmer(bks1)
Anova(bks1)
lsmeans(bks1, pairwise~Species, adjust="tukey")

bks2 <- glmer.nb(NatHerb_Rich ~ Species + (1|Block) , data = spray)
summary(bks2)
plot(bks2)
dispersion_glmer(bks2)
Anova(bks2)
lsmeans(bks2, pairwise~Species, adjust="tukey")


#contrasts
summary(glht(bks1, linfct = mcp(Species = "Tukey")))


######Bark; response: Exotic Abundance
hist(bark$SDI)

bkiv0 <- lmer(SDI ~ Species + (1|Block) , data = bark)
summary(bkiv0)
plot(bkiv0)
Anova(bkiv0)
AIC(bkiv0)
lsmeans(bkiv0, pairwise~Species, adjust="tukey")

bkiv1 <- glmer(SDI ~ Species + (1|Block) , family=Gamma, data = bark)
summary(bkiv1)
plot(bkiv1)
dispersion_glmer(bkiv1)
Anova(bkiv1)
AIC(bkiv1)
anova(bkiv1,bkiv0)
lsmeans(bkiv1, pairwise~Species, adjust="tukey")

bkiv2 <- glmer.nb(Invasive ~ Species +  (1|Block),family = inverse.gaussian, data = bark)
summary(bkiv2)
plot(bkiv2)
dispersion_glmer(bkiv2)
Anova(bkiv2)

anova(bkiv2,bkiv1)

bkiv3 <- glmer(SDI ~ Species +  (1|Block), family=inverse.gaussian, data = bark)
summary(bkiv3)
plot(bkiv3)
dispersion_glmer(bkiv3)
Anova(bkiv3)

anova(bkiv3, bkiv1)

#contrasts
summary(glht(bkiv1, linfct = mcp(Species = "Tukey")))
summary(glht(bkiv2, linfct = mcp(Period = "Tukey")))

K5 <- glht(bkiv2, linfct = mcp(Species = "Tukey"))$linfct
K6 <- glht(bkiv2, linfct = mcp(Period = "Tukey"))$linfct

summary(glht(bkiv2, linfct = rbind(K3, K4)))

#####
#New day, Just trying Species effects
#####Bark; response: Total abundance
hist(spray$SDI)
spray$native <- (spray$Tot_abun - spray$Invasive)

ta1 <- lmer(SDI ~ Species + (1|Block), data = spray)
summary(ta1)
plot(ta1)
dispersion_glmer(ta1)
Anova(ta1)

qqnorm(resid(ta1))
qqline(resid(ta1))

ta2 <- glmer(SDI ~ Species +  (1|Block),family = Gamma, data = spray)
summary(ta2)
plot(ta2)
dispersion_glmer(ta2)
Anova(ta2)

anova(ta2, ta1)

lsmeans(ta2, pairwise~Species, adjust="tukey")

#contrast
summary(glht(ta2, linfct = mcp(Species = "Tukey")))



###################################Herbivory Data 6.4.2020#################
#Enter data
herb <- read.csv("herbdata.csv", header=TRUE)
herb

#Set Period, Block, Tree and year as factors vs. integers
herb$Plot <- as.factor(herb$Plot)

hist(herb$logit2)


lm1 <- lmer(logit~Species + (1|Plot), data = herb)
summary(lm1)
plot(lm1)

lm2 <- lmer(logit2~Species+Month +  (1|Plot/Tree), data = herb)
summary(lm2)
plot(lm2)
AIC(lm2)
qqnorm(resid(lm2))
qqline(resid(lm2))

lm3 <- lmer(logit2~Species* Month+ (1|Plot/Tree), data = herb)
summary(lm3)
plot(lm3)
Anova(lm3)
AIC(lm3)
qqnorm(resid(lm3))
qqline(resid(lm3))

relgrad <- with(lm3@optinfo$derivs,solve(Hessian,gradient))
max(abs(relgrad))

anova(lm3,lm2)

#contrasts
summary(glht(lm3, linfct = mcp(Species = "Tukey")))
summary(glht(lm3, linfct = mcp(Month = "Tukey")))

K5 <- glht(lm3, linfct = mcp(Species = "Tukey"))$linfct
K6 <- glht(lm3, linfct = mcp(Month = "Tukey"))$linfct

summary(glht(lm3, linfct = rbind(K5, K6)))

herb$ms <- with(herb, interaction(Month, Species))
cell <- lmer(logit2~ ms + (1|Plot/Tree), data = herb)
summary(glht(cell, linfct = rbind(K5,K6)))

install.packages("lsmeans")
library(lsmeans)
lsmeans(lm3, pairwise~Month*Species, adjust="tukey")

citation("lme4")
citation("emmeans")
citation("car")
citation("blmeco")
