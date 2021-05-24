library(tidyr)
library(dplyr)
library(mvabund)
library(ggplot2)
library(Rmisc)
library(patchwork)
library(MuMIn)
library(vegan)

#FISH ANALYSES
#pivot table the data so columns are species


All_data <- read_csv("C:/Users/vr371/Dropbox/Papers/Vanuatu depth range/All data fishes cleaned.csv")

Fishdat <-  All_data %>% 
  pivot_wider(names_from = SpeciesFull, values_from = MaxN, values_fill = 0)


Fishdat$Depth <- as.factor(Fishdat$Depth)


#start mvabund analyses

Fish <- Fishdat[,4:108]

Fish.sq <- sqrt(Fish)

Fishstuff <- mvabund(Fish.sq)


par(mar=c(2,10,2,2)) # adjusts the margins
boxplot(Fishstuff[,5:11],horizontal = TRUE,las=2, main="Abundance")

#check variances

meanvar.plot(Fishstuff)

#model
#try different transformations and check residuals


mod1 <- manyglm(Fishstuff ~ Fishdat$DayNight*Fishdat$Depth, family = "negative_binomial")
mod2 <- manyglm(Fishstuff ~ Fishdat$DayNight*Fishdat$Depth, family = "poisson")

plot(mod1)
plot(mod2)


#less fan shape for poisson distribution, so use that

anova(mod2)


#Pairwise analysis

anova(mod2, p.uni = "adjusted")


#NMDS plot for fish communities

#remove species/rows that are not informative

Fishdat <- read_csv("Fishdat.csv")

test <- metaMDS(Fishdat[,5:22], k = 2, trymax = 1000)

fishnmds <- metaMDS(Fishdat[,5:22], k = 2, trymax = 500, previous.best = test)

stressplot(fishnmds)

fishnmds$stress

data.scores <- as.data.frame(scores(fishnmds))  #Using the scores function from vegan to extract the site scores and convert to a data.frame
data.scores$DayNight <- Fishdat$DayNight  #  add the grp variable created earlier
data.scores$Depth <- Fishdat$Depth

ggplot(data = data.scores, aes(x = NMDS1, y = NMDS2)) + geom_point(shape = 16, aes(size = Depth, colour = DayNight)) + theme_classic() + 
  scale_colour_grey() + theme(axis.text.y = element_blank(), axis.text.x = element_blank())

#Cumulative abundance & Species Richness

Fishdat$Cumulative <- rowSums(Fishdat[,4:108])

Fishdat$Richness <- rowSums(Fishdat[,4:108] > 0)

Fishdat$DayNumber <- as.factor(Fishdat$DayNumber)

meanCumulative <- summarySE(data = Fishdat, measurevar = "Cumulative", groupvars = c("Depth", "DayNight"))

meanRichness <- summarySE(data = Fishdat, measurevar = "Richness", groupvars = c("Depth", "DayNight"))

meanCumulative$Depth <- as.factor(meanCumulative$Depth)


#models for fish stuff

maxnmod <- glm(Cumulative ~ Depth*DayNight, family = Gamma(link = "inverse"), data = Fishdat)

summary(maxnmod)

richmod <- glm(Richness ~ Depth*DayNight, family = Gamma(link = "inverse"), data = Fishdat)

summary(richmod)



#Graphs

a <- ggplot(data = Fishdat, aes(y = Cumulative, x = Depth)) + geom_point(aes(colour = DayNight)) + theme_classic() + geom_smooth(aes(colour = DayNight)) +
  ylab("Cumulative MaxN") + xlab("Depth (m)") + scale_colour_grey()

b <- ggplot(data = Fishdat, aes(y = Richness, x = Depth)) + geom_point(aes(colour = DayNight)) + theme_classic() + geom_smooth(aes(colour = DayNight)) +
  ylab("Species Richness") + xlab("Depth (m)") + scale_colour_grey()

a/b

#substrate data

Substrates <- read_csv("C:/Users/vr371/Dropbox/Papers/Vanuatu depth range/Substrates.csv")

Substrates$Substrate_type <- as.factor(Substrates$Substrate_type)

substratemod <- glm(Substrate_type ~ Depth, data = Substrates, family = binomial)

Substrates$Depth <- as.factor(Substrates$Depth)

ggplot(data = Substrates, aes(x = Depth)) + geom_bar(position = "fill", aes(fill = Substrate_type), colour = "black") + theme_classic() + xlab("Depth (m)") +
  scale_fill_viridis_d() + ylab("Proportion of measures") + labs(fill = "Substrate class")



#Benthic organism data

Benthic_organisms <- read_csv("C:/Users/vr371/Dropbox/Papers/Vanuatu depth range/Benthic organisms.csv")


#assemblage analysis for benthic stuff


Benthicstuff <- mvabund(Benthic_organisms[,3:7])


mod3 <- manyglm(Benthicstuff ~ Benthic_organisms$DayNight*Benthic_organisms$Depth, family = "negative_binomial")
mod4 <- manyglm(Benthicstuff ~ Benthic_organisms$DayNight*Benthic_organisms$Depth, family = "poisson")
plot(mod3)
plot(mod4)

#negative binomial better

anova(mod3)

#visualize data

Benthicdat <-  Benthic_organisms %>% 
  pivot_longer(cols = SoftCoral:Fan, names_to = "Organisms", values_to = "Count")

Benthicmeans <- summarySE(data = Benthicdat, measurevar = "Count", groupvars = c("Depth", "Organisms"))

Benthicmeans$Depth <- as.factor(Benthicmeans$Depth)

ggplot(data = Benthicmeans, aes(x = Depth, y = Count)) + geom_col() + theme_classic() + xlab("Depth (m)") + 
  facet_wrap(DayNight~ Organisms) + ylab("Mean count per transect") + geom_errorbar(aes(ymax = Count + se, ymin = Count - se), width = 0.5)


