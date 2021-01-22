library(vegan)
library(lme4)
library(car)
library(performance)

# Fitting Maximum‐likelihood population‐effects (MLPE) mixed models in lme4 using population identity comparisons as a random effect
dat <- read.table("genotype_matrix.txt")
SNP <- as.data.frame(t(dat))

env = read.table("environmental.txt", header=TRUE)
EnvVars <- as.data.frame(scale(env[, 5:9], scale=TRUE, center=TRUE))
landscape <- as.data.frame(scale(env[, 13:16], scale=TRUE, center=TRUE))
space <- as.data.frame(scale(env[,3:4], scale=TRUE, center=TRUE))


env_dist_matrix <-  dist(EnvVars, upper=TRUE, method = "euclidean")
env_dist <- as.matrix(env_dist_matrix, upper=TRUE, diagonal = TRUE)
env_dist <- ResistanceGA:::lower(env_dist)

gen_dist <- as.matrix(dist(SNP, upper = TRUE))
gen_dist <- ResistanceGA:::lower(gen_dist)

geo_dist_matrix <-  dist(space, upper=TRUE, method = "euclidean")
geo_dist <- as.matrix(geo_dist_matrix, upper=TRUE, diagonal = TRUE)
geo_dist <- ResistanceGA:::lower(geo_dist)
                      
land_dist_matrix <-  dist(landscape, upper=TRUE, method = "euclidean")
land_dist <- as.matrix(land_dist_matrix, upper=TRUE, diagonal = TRUE)
land_dist <- ResistanceGA:::lower(land_dist)
                      
pop <- read.table("pop.txt", header = TRUE)
                      
# null model
null <- lmer(gen_dist ~ 1  + (1 | pop$pop1), REML = T)
summary(null)
Anova(null)
r2(null)
                      
# Model1 - gen_dist ~ geo_dist (IBD)
IBD <- lmer(gen_dist ~ geo_dist  + (1 | pop$pop1), REML = T)
summary(IBD)
Anova(IBD)
r2(IBD)
                      
# Model2- gen_dist ~ env_dist (IBE)
IBE <- lmer(gen_dist ~ env_dist  + (1 | pop$pop1),REML = T)
summary(IBE)
Anova(IBE)
r2(IBE)
                      
# Model3 - gen_dist ~ land_dist (~IBR)
IBR <- lmer(gen_dist ~ land_dist  + (1 | pop$pop1), REML = T)
summary(IBR)
Anova(IBR)
r2(IBR)
                      
# Model4 - gen_dist ~ geo_dist + env_dist (IBD + IBE)
IBD_IBE <- lmer(gen_dist ~ geo_dist*env_dist  + (1 | pop$pop1), REML = T)
summary(IBD_IBE)
Anova(IBD_IBE)
r2(IBD_IBE)
                      
# Model5 - gen_dist ~ geo_dist + land_dist (IBD + ~IBR)
IBD_IBR <- lmer(gen_dist ~ geo_dist*land_dist  + (1 | pop$pop1), REML = T)
summary(IBD_IBR)
Anova(IBD_IBR)
r2(IBD_IBR)

# Model6 - gen_dist ~ env_dist + land_dist (IBE + ~IBR)
IBE_IBR <- lmer(gen_dist ~ env_dist*land_dist  + (1 | pop$pop1), REML = T)
summary(IBE_IBR)
Anova(IBE_IBR)
r2(IBE_IBR)

AIC(null, IBD, IBE, IBR, IBE_IBR, IBD_IBE, IBD_IBR)

anova(IBD_IBE, IBD_IBR)
