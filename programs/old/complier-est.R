#!/usr/bin/R
## Senan Hogan-Hennessy, 28 March 2023
## Simulation for estimating compliers.
print(Sys.time())
set.seed(47)
# Functions for data manipulation and visualisation
library(tidyverse)
# Generalised Random Forests https://grf-labs.github.io/grf/
library(grf)
# Functions for fast linear models with IV + FEs
library(lfe)
# Functions for hypothesis testing
library(car)
library(lmtest)
# TeX tables
library(stargazer)
# Define number of digits in tables and graphs
digits.no <- 3
# Size for figures
fig.width <- 10
fig.height <- fig.width * 0.85


################################################################################
## Simulate data.

# Sample size
sample.size <- 10^3
# Number of dimensions for the covariates
covariate.dim <- 5
# Generate matrix of covariates, along which treatment intensity varies
min_X <- 0
max_X <- 1 / covariate.dim
X <- runif(sample.size * covariate.dim, min = min_X, max = max_X)
X <- matrix(X, sample.size, covariate.dim)
# Specifiy the treatment heterogeneity.  Start with everyone being non-compliers
firststage.designation <- replicate(sample.size, "never-taker", "vector")
pi_X <- replicate(sample.size, 0, "vector")
# some are never-takers pi(x) = 0, compliers pi(X) = 1, always-takers pi(X) = 0
complier.threshold <- 5 / 12
alwaystaker.threshold <- 7 / 12
pi_X[rowSums(X) < complier.threshold] <- 0
pi_X[rowSums(X) >= complier.threshold] <- 1
firststage.designation[rowSums(X) > complier.threshold] <-
    replicate(sum(rowSums(X) > complier.threshold), "complier", "vector")
pi_X[rowSums(X) >= alwaystaker.threshold] <- 0
firststage.designation[rowSums(X) > alwaystaker.threshold] <-
    replicate(sum(rowSums(X) > alwaystaker.threshold), "always-taker", "vector")
# Plot the first stage heterogeneity
plot(rowSums(X), pi_X)
hist(pi_X)
# Generate the list of a randomly assigned instrument.
probZ <- 0.5
Z <- rbinom(sample.size, 1, probZ)
# Generate the list of a binary treatment, take-up heterogeneity by pi(X)
# Treatment propensity, determined by compliers score pi(X) and instrument Z
a <- ifelse(rowSums(X) >= alwaystaker.threshold, 1, 0)
D <- a + pi_X * Z
# Show difference in first-stage LATE vs ATE.
lm(D[firststage.designation == "complier"] ~
    1 + Z[firststage.designation == "complier"]) %>%
    summary()
lm(D ~ 1 + Z) %>% summary()
# Generate the outcome, with het treatment effect tau(x)
tau_X <- replicate(sample.size, 0, "vector")
tau_X[rowSums(X) <= complier.threshold] <-
    rnorm(sum(rowSums(X) <= complier.threshold), 1, 0.1)
tau_X[rowSums(X) > complier.threshold] <-
    rnorm(sum(rowSums(X) > complier.threshold), 2, 0.1)
tau_X[rowSums(X) > alwaystaker.threshold] <-
    rnorm(sum(rowSums(X) > alwaystaker.threshold), 3, 0.1)
Y <- tau_X * D + rnorm(sample.size, 0, 0.2)
# Store the real LATE & ATE, for the first stage.
firststage_ATE <- mean(pi_X)
firststage_LATE <- mean(pi_X[firststage.designation == "complier"])
print(c(firststage_ATE, firststage_LATE))
# Store the real LATE & ATE, for the second stage.
secondstage_ATE <- mean(tau_X)
secondstage_LATE <- mean(tau_X[firststage.designation == "complier"])
# Put these data to a coherent dataframe.
combined.data <- tibble(
    # Observed data (minus ambiguous X column length)
    Z = Z, D = D, Y = Y,
    # Unobserved.
    latent_X = rowSums(X),
    complier = firststage.designation,
    pi_X = pi_X,
    tau_X = tau_X)


# Does controlling for a mediator have an effect?
# Suppose we have an idealised experiment for Z, but a mediating causal effect:
# Z -> T -> Y, and also direct effect Z -> Y
# Does Z -> Y become identified if we control for T, or is this collider bias?
#! ANswer: Z-> Y becomes identified in special cases, unsure if general.

lm(Y ~ 1 + Z) %>% summary()
mean(tau_X * pi_X)
secondstage_ATE
lm(Y ~ 1 + predict(lm(D ~ 1 + Z))) %>% summary()
secondstage_LATE




################################################################################
## Estimate complier score in the first stage, via MLE ICSW

# Estimate the complier score with Arnow Carnegie (2013) ICSW,
# which is a special case of Abadie (2003) kappa weighting.
library(icsw)
cscoreout <- compliance.score(D = D, Z = Z, W = X, genoud = TRUE)
# Extract vector of estimated compliance scores
combined.data$icsw_complier_score <- cscoreout$C.score
# SHow the estimates of complier score in a figure.
icsw_complier.plot <- combined.data %>%
    ggplot(aes(x = icsw_complier_score)) +
    geom_histogram(aes(fill = complier), alpha = 0.75) +
    labs(fill = "Oracle Designation:") +
    theme_bw() +
    theme(legend.position = "top")
icsw_complier.plot
# SHow the 2SLS efficency gain by complier weighting.
ivreg::ivreg(Y ~ D | Z, data = combined.data) %>% summary()
ivreg::ivreg(Y ~ D | Z, weights = combined.data$icsw_complier_score,
    data = combined.data) %>% summary()


################################################################################
## Estimate Heterogenous TE in the first stage, via GRF

# Estimate the first stage by GRF causal forest
count.trees <- 2000
firststage.forest <- causal_forest(X, D, Z,
    tune.parameters = "all",
    num.trees = count.trees,
    tune.num.trees = count.trees)
firststage.pred <- predict(firststage.forest, X, estimate.variance = TRUE)
combined.data$firststage_pred <- firststage.pred$predictions
# Show the covariate determinant and first-stage estimates, and real first-stage
combined.data %>%
    ggplot(aes(x = latent_X, colour = complier)) +
    geom_point(aes(y = pi_X), alpha = 0.75) +
    geom_point(aes(y = firststage_pred), alpha = 0.75) +
    labs(colour = "Oracle Designation:") +
    theme_bw() +
    theme(legend.position = "top")


################################################################################
## Estimate complier score, by GRF for Abadie (2003) complier score.

count.trees <- 2000
# Abadie (2003) proves Pr(D(1) > D(0) | X) = E[D | X, Z=1] - E[D | X, Z=0]
# With similar results for ATs + NTs, by the same result.
# Estimate Pr(AT | X) = E[D | X, Z = 0] by honest GRF
alwaystaker.forest <- regression_forest(X[Z==0, ], D[Z==0],
    tune.parameters = "all",
    num.trees = count.trees)
combined.data$alwaystaker_score <- predict(alwaystaker.forest, X)$predictions
# Estimate Pr(NT | X) = 1 - E[D | X, Z = 1] by honest GRF
nevertaker.forest <- regression_forest(X[Z==1, ], D[Z==1],
    tune.parameters = "all",
    num.trees = count.trees)
combined.data$nevertaker_score <- 1 - predict(nevertaker.forest, X)$predictions
# Combine for the complier score.
combined.data$complier_score <- (1 - combined.data$nevertaker_score -
    combined.data$alwaystaker_score)
# Show the covariate determinant and first-stage estimates, and real first-stage
combined.data %>%
    ggplot(aes(x = latent_X, colour = complier)) +
    geom_point(aes(y = pi_X), alpha = 0.75) +
    geom_point(aes(y = complier_score), alpha = 0.75) +
    labs(colour = "Oracle Designation:") +
    theme_bw() +
    theme(legend.position = "top")

# SHow the 2SLS efficency gain by complier weighting.
ivreg::ivreg(Y ~ D | Z, data = combined.data) %>% summary()
ivreg::ivreg(Y ~ D | Z, weights = complier_score,
    data = combined.data) %>% summary()

# This method works much better, and so can be used to profile compliers much better.
mean(combined.data[combined.data$complier == "complier", ]$abadie_complier_score)
mean(combined.data[combined.data$complier == "never-taker", ]$abadie_complier_score)
mean(combined.data[combined.data$complier == "always-taker", ]$abadie_complier_score)

## Test the exclusion restriction among definite non-compliers
combined.data %>%
    #filter(complier_score < 0.5) %>%
    ggplot(aes(x = complier_score)) +
    geom_histogram()
combined.data %>%
    filter(complier_score < 0.5) %>%
    lm(Y ~ 1 + Z + D, data = .) %>%
    summary()
#! THE EXCLUSION RESTRICTION TEST WORKS REALLY, REALLY WELL.


################################################################################
## Latent index function learning

# Estimate the implied latent function by prob forest.
latent.forest <-
    probability_forest(X, factor(D), num.trees = 10000)
combined.data$latent_predictions <-
    predict(latent.forest, X)$predictions[, 2]
# Show the correspondence.
combined.data %>%
    ggplot(aes(x = latent_X, y = latent_predictions)) +
    geom_point() +
    stat_smooth() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = "red")
# And plot the treatment take-up along estimated latemt functionm
hist(combined.data$latent_predictions)

# Show how E[D | X] varies with real V(Z; X)
combined.data %>%
    ggplot(aes(x = latent_X, y = D, colour = Z)) +
    geom_point() +
    geom_smooth(span = 1 / 10, colour = Z)

#! Test the score estimation
secondstage.forest <- instrumental_forest(X, Y, D, Z,
    tune.parameters = "all",
    num.trees = count.trees)
combined.data$scores <- get_scores(secondstage.forest)
mean(combined.data$scores)
mean(tau_X[firststage.designation == "complier"])
mean(tau_X[firststage.designation == "never-taker"])
mean(tau_X[firststage.designation == "always-taker"])


#! Additionally, how do I learn the threshold function latent_X, i.e. V(Z; X)
#! WHich here is equal to \sum X across the columns.
#! I can, in principle, learn this function as a way of decribing + explaining
#! programme participation

# Estimator road-map
# 1. Causal forest to detect heterogeneity in first stage, along observed X.
#    Proceed with a decision rule for detecting firststage-het along observed X.
# 2. Apply the Abadie complier score method to profile compliers + NTs + ATs
# 3. Validate exclusion + monotonicity using the complier score
#   -> If these are violated, provide estimators that correct for this
#   (Cousins notes the LATE estimand needs refining if there are defiers.)
# 4. Apply the Abadie+ (2022) LATE efficiency gain with forest chosen X values.
# 5. Reweight the LATE to ATE (following Angrist Fer 2013, Aronow+ 2013)

#TODO: write closed form expressions for following, using Abadie (2003) approach
#TODO: NT score, Pr(D(1) = D(0) = 0 | X)
#TODO: AT score, Pr(D(1) = D(0) = 1 | X)
#TODO: defier score, Pr(D(1) < D(0) | X)
#! -> defier score is hard to identify, though it may be possible with policy
#!    learning https://grf-labs.github.io/grf/articles/policy_learning.html
#!    i.e. find the score for a negative local-effect

#TODO: Think of an empirical example which will be perfect to tailor this for,
#TODO: likely in a higher education setting, with available covariate info.
#TODO: Perhaps available within UK data.
