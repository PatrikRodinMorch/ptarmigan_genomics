library(readr)
library(vroom)
library(vegan)

# Inspired by Forester et al. 2018. Mol Ecol.

# In this analysis we are aiming to find association between allele frequencies in the different populations
# and the environmental variables that individuals in those populations experience.
# For that we will use a redundancy analysis (RDA) which can selection imposed by multivariate environmental gradients
# The RDA in essence consists of two tests 1) a multiple regression with a matrix of allele counts as response variables and 4-6 environmental variables as predictors
# 2) a PCA is then run on the fitted values producing linear combination of the predictors.

# Read in the environmental data for both rock ptarmigan and willow ptarmigan.
# These data have been cleaned up and variables showing high correlation have been dropped using:
# data_prep_envass.R

muta_env <- read_csv("muta_env.csv")
lago_env <- read_csv("lago_env.csv")

# Standardize the environmental variables to account for the difference in scale of the measurement of each variable
muta_env[, 3:10] <- scale(muta_env[, 3:10], scale = TRUE, center = FALSE)
lago_env[, 3:8] <- scale(lago_env[, 3:8], scale = TRUE, center = FALSE)


# Read in the genome wide data to be used as response variables. 
# The input is in plink raw format and coded 0/1/2 for homozygots for each allele and heterozygots
# Since these files are large (~ 1.5 million columns and 61 rows) we need to use something different than base read.cvs or read_csv from the readr package.
# vroom is a great package for reading in larger files.

muta_gen <- vroom("muta_filtered_envass_noheader.raw", col_names = FALSE, delim = " ", num_threads = 12)
lago_gen <- vroom("lago_filtered_envass_noheader.raw", col_names = FALSE, delim = " ", num_threads = 12)

# The RDA models
muta_rda <- rda(muta_gen[, 7:length(muta_gen)] ~ elevation + bio5 + bio6 + bio8 + bio15, data=muta_env)
muta_rda
lago_rda <- rda(7:length(lago_gen) ~ elevation + bio5 + bio6 + bio7 + bio8 + bio15, data=lago_env)
lago_rda

# The variance explained by the environmental predictors (constrained axes)
muta_rda_r2 <- RsquareAdj(muta_rda)
muta_rda_r2
muta_rda_eigenvals <- summary(eigenvals(muta_rda, model = "constrained")) # looking closer at each constrained axis
screeplot(muta_rda) # Visualizing the eigenval summary output

lago_rda_r2 <- RsquareAdj(lago_rda)
lago_rda_eigenvals <- summary(eigenvals(lago_rda, model = "constrained"))
screeplot(lago_rda)

# Test partial RDA which excludes the influence of some predictor. In this case we want to correct for spatial autocorrelation in allele frequencies.
# We might have to control for spatial autocorrelation (population structure) in some other way
muta_prda <- rda(muta_gen[, 7:length(muta_gen)] ~ elevation + bio5 + bio6 + bio7 + bio8 + bio15 + Condition(Long + Lat), data=muta_env)
muta_prda
lago_prda <- rda(lago_gen[, 7:length(lago_gen)] ~ elevation + bio5 + bio6 + bio7 + bio8 + bio15 + Condition(Long + Lat), data=lago_env)
lago_prda

# The variance explained by the environmental predictors (constrained axes)
muta_prda_r2 <- RsquareAdj(muta_prda)
muta_prda_r2
muta_prda_eigenvals <- summary(eigenvals(muta_prda, model = "constrained")) # looking closer at each constrained axis
screeplot(muta_prda) # Visualizing the eigenval summary output

lago_prda_r2 <- RsquareAdj(lago_prda)
lago_prda_eigenvals <- summary(eigenvals(lago_prda, model = "constrained"))
screeplot(lago_prda)

# Testing the significance of the full model using the ANOVA like permutation test anova.cca
sign_full_model_muta <- anova.cca(muta_rda, parallel = getOption("mc.cores"), permutation = 999)
sign_full_model_lago <- anova.cca(lago_rda, parallel = getOption("mc.cores"), permutation = 999)

# We can also asses the significance of each axis, somethign that is extremely computationally intensive and takes a long time
sign_axis_muta <- anova.cca(muta_rda, by = "axis", parallel = getOption("mc.cores"))
sign_axis_lago <- anova.cca(lago_rda, by = "axis", parallel = getOption("mc.cores"))

# Plotting the results of the full model with symmetrical scaling
pop_col <- c("#ff7f00", "#fff200", "#2bff00", "#00ffee", "#3711f5", "#f511ed", "#f51115")

# RDA axis 1 and 2
plot(muta_rda, type="n", scaling=3)
points(muta_rda, display="species", pch=20, cex=0.7, col="gray32", scaling=3)           # the SNPs
points(muta_rda, display="sites", pch=21, cex=1.3, col="gray32", scaling=3, bg=pop_col[as.factor(muta_env[["Pop"]])]) # the birds
text(muta_rda, scaling=3, display="bp", col="#0868ac", cex=1, choices=c(1,2))                          # the environmental variables
legend("bottomright", legend=levels(as.factor(muta_env[["Pop"]])), bty="n", col="gray32", pch=21, cex=1, pt.bg=pop_col)


library(ggord)
ggord(muta_rda, as.factor(muta_env[["Pop"]])) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# Now we use the method of Forester et al. 2018 Mol.Ecol to identify candidates involved in adaptation
# The method extract the loadings for each RDA axis then identify the genetic variants (SNPs) whose loadings falls within either tail of the distribution of loading values
loads.muta_rda <- scores(muta_rda, choices = 1:3, display = "species")

ggplot(aes(RDA1), data =data.frame(loads.muta_rda)) + geom_histogram( binwidth = 0.005) + ggtitle("RDA1 loadings")
ggplot(aes(RDA2), data =data.frame(loads.muta_rda)) + geom_histogram( binwidth = 0.005) + ggtitle("RDA2 loadings")
ggplot(aes(RDA3), data =data.frame(loads.muta_rda)) + geom_histogram( binwidth = 0.005) + ggtitle("RDA3 loadings")


# This function loadings in a vector and the number of SD from the mean where we want to set our cutoff to identify candidate SNPs
candidates <- function(loads, z) {
  lims <- mean(loads) + c(-1, 1) * z * sd(loads)
  loads[loads < lims[1] | loads > lims[2]]
} 

# Save the candidates identified to be at least |3| standard deviations from the mean, for each RDA axis in separate vectors then count the total number
cand1 <- candidates(loads.muta_rda[ ,"RDA1"], 3)
cand2 <- candidates(loads.muta_rda[ ,"RDA2"], 3)
cand3 <- candidates(loads.muta_rda[ ,"RDA3"], 3)

ncand <- length(cand1) + length(cand2) + length(cand3)
ncand
