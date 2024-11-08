#!/usr/bin/R
## Senan Hogan-Hennessy, 06 February 2023
## Simulation for hetergogenous IV estimation.
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
sample.size <- 10^4
# Number of dimensions for the covariates
covariate.dim <- 10
# Generate matrix of covariates, along which treatment intensity varies
min_X <- 0
max_X <- 1 / covariate.dim
X <- runif(sample.size * covariate.dim, min = min_X, max = max_X)
X <- matrix(X, sample.size, covariate.dim)
# Specifiy the treatment heterogeneity.  Start with everyone being non-compliers
firststage.designation <- replicate(sample.size, "non-complier", "vector")
pi_X <- replicate(sample.size, 0, "vector")
# some are non-compliers pi(x) = 0, some are compliers pi(X) = 1
complier.threshold <- 0.5
pi_X[rowSums(X) > complier.threshold] <-
    rnorm(sum(rowSums(X) <= complier.threshold), 0.25, 0.1)
pi_X[rowSums(X) > complier.threshold] <-
    rnorm(sum(rowSums(X) > complier.threshold), 0.75, 0.1)
firststage.designation[rowSums(X) > complier.threshold] <-
    replicate(sum(rowSums(X) > complier.threshold), "complier", "vector")
#! Plot the first stage heterogeneity
plot(rowSums(X), pi_X)
hist(pi_X)
# Generate the list of a random instrument.
probZ <- 0.5
Z <- rbinom(sample.size, 1, probZ)
# Generate the list of a binary treatment, take-up heterogeneity by pi(X)
# Treatment propensity, determined by compliers score pi(X) and instrument Z
D <- rbinom(sample.size, 1, pi_X * Z)
# Generate the outcome, with treatment effect tau(x)
tau_X <- replicate(sample.size, 0, "vector")
tau_X[rowSums(X) <= complier.threshold] <-
    rnorm(sum(rowSums(X) <= complier.threshold), 0.25, 0.1)
tau_X[rowSums(X) > complier.threshold] <-
    rnorm(sum(rowSums(X) > complier.threshold), 0.75, 0.1)
Y <- tau_X * D + rnorm(sample.size, 0, 0.1)
# Store the real LATE & ATE, for the first stage.
firststage_ATE <- mean(pi_X)
firststage_LATE <- mean(pi_X[firststage.designation == "complier"])
#! SHowfirst stage ATE + LATE
print(c(firststage_ATE, firststage_LATE))
# Store the real LATE & ATE, for the second stage.
secondstage_ATE <- mean(tau_X)
secondstage_LATE <- mean(tau_X[firststage.designation == "complier"])
# Put these data to a coherent dataframe.
combined.data <- data.frame(
    # Observed data
    Z = Z, D = D, Y = Y,
    # Unobserved.
    "complier?" = firststage.designation,
    pi_X = pi_X, tau_X = tau_X)
#! NOTE: only never-takers exist in this DGP
#TODO: adjust to include always-takers.


################################################################################
## Estimate the heterogenous TE in the first stage by GRF.

# Estimate the first stage by GRF causal forest
count.trees <- 50
firststage.forest <- causal_forest(X, D, Z, num.trees = count.trees)
firststage.pred <- predict(firststage.forest, X, estimate.variance = TRUE)
firststage.forest

# Estimate complier score within each tree of the pre-built RF
matrix.complier.tstats <- matrix(NA, sample.size, count.trees)
matrix.complier.pvals <- matrix(NA, sample.size, count.trees)
for (count.tree in 1:count.trees){
    print(count.tree)
    est.complier.tstats <- rep(NA, sample.size)
    est.complier.pvals <- rep(NA, sample.size)
    tree <- get_tree(firststage.forest, count.tree)
    nodes <- get_leaf_node(tree, X)
    for (node in unique(nodes)){
        tryCatch({
            complier.test <- lm(D[nodes == node] ~ Z[nodes == node])
            complier.tstat <- summary(complier.test
                )[["coefficients"]]["Z[nodes == node]", "t value"]
            complier.pval <- summary(complier.test
                )[["coefficients"]]["Z[nodes == node]", "Pr(>|t|)"]
            est.complier.tstats[nodes == node] <-
                rep(complier.tstat, sum(nodes == node))
            est.complier.pvals[nodes == node] <-
                rep(complier.pval, sum(nodes == node))
        }, error = function(e) {
            est.complier.tstats[nodes == node] <-
                rep(NA, sum(nodes == node))
            est.complier.pvals[nodes == node] <-
                rep(NA, sum(nodes == node))
        })
    }
    # Save the results
    matrix.complier.tstats[, count.tree] <- est.complier.tstats
    matrix.complier.pvals[, count.tree] <- est.complier.pvals
}
# Average across the trees
complier.tstats <- rowMeans(matrix.complier.tstats, na.rm = TRUE)
complier.pvals <- rowMeans(matrix.complier.pvals, na.rm = TRUE)

hist(complier.pvals)
#! PROBLEM -> COMPLIER SCORES IS STILL OFF
#! I WANT A HUMP AT AROUND 0 FOR NON-COMPLIERS
#! AND A HUMP AT 1 FOR COMPLIERS.


# Try the IV Forest
secondstage.forest <- instrumental_forest(X, Y, D, Z,
    num.trees = count.trees)
average_treatment_effect(
    secondstage.forest, target.sample = "all")
# Try the IV Forest with weights for complier score.
secondstage.forest <- instrumental_forest(X, Y, D, Z,
    sample.weights = 1 - complier.pvals,
    num.trees = count.trees)
average_treatment_effect(
    secondstage.forest, target.sample = "all")
#! SHowfirst stage ATE + LATE
print(c(secondstage_ATE, secondstage_LATE))

#! Weighting scheme does not work -> Try the approach of rejecting the null within each tree instead.




instrumental_forest(X, Y, W, Z)

#! THIS APPROACH WORKS!  IT WORKS, BATMAN!
average_treatment_effect(
    firststage.forest, target.sample = "all",
    debiasing.weights = rep(1, sample.size))
average_treatment_effect(
    firststage.forest, target.sample = "all",
    debiasing.weights = rep(1, sample.size),
    subset = (complier.pvals == 1))
average_treatment_effect(
    firststage.forest, target.sample = "all",
    debiasing.weights = (1 / complier.tstats))


# Estimate the conditional average treatment effect on the full sample (CATE).
est_firststage_ATE <- average_treatment_effect(
    firststage.forest, target.sample = "all")
# Estimate the conditional average treatment effect on compliers
est_firststage_LATE <- average_treatment_effect(
    firststage.forest, target.sample = "all",
    subset = (firststage.designation == "complier"))
# Compare estimates to real (with oracle knowledge of pi(X)) -> they are correct
print(c("real:", firststage_ATE, "estimate:", est_firststage_ATE))
print(c("real:", firststage_LATE, "estimate:", est_firststage_LATE))
# Treatment effect among non-compliers.
average_treatment_effect(
    firststage.forest, target.sample = "all",
    subset = (firststage.designation == "non-complier"))

#! Problem: the above LATE estimate uses oracle knowledge of complier status.
# Try to get CATE estimates of where \hat \pi(X) > 0
firststage.pred <- predict(firststage.forest, X, estimate.variance = TRUE)
hist(firststage.pred$predictions)
hist(firststage.pred$predictions[firststage.designation == "complier"])
table(pi_X[firststage.designation == "complier"])
hist(firststage.pred$predictions[firststage.designation == "non-complier"])
hist(pi_X)
plot(rowSums(X), pi_X)
# Plot the individual estimated TE.
plot(rowSums(X), firststage.pred$predictions)



#! PACKAGE THAT DOES THIS ESTIMATION:
# https://rdrr.io/cran/icsw/man/compliance.score.html
# install.packages("icsw")
library(icsw)
cscoreout <- compliance.score(D = D, Z = Z, W = X, genoud = TRUE)
# Extract vector of estimated compliance scores
complier.score <- cscoreout$C.score
complier.prediction <- as.integer(complier.score > 0.5)
table(complier.prediction, firststage.designation)
# Try with model miss-specification -> complier classification does not work.
cscoreout <- compliance.score(D = D, Z = Z,
    W = X[, 1 : (covariate.dim - 3)], genoud = TRUE)
complier.score <- cscoreout$C.score
noncomplier.prediction <- as.integer(complier.score <= 0.5)
table(noncomplier.prediction, firststage.designation)

#! ARNOW 2013 WORKS!  (See Kennedy+ 2020 for conditions on testable cond on \pi(X)
#! That allows for classification of compliers).
#! I wonder if this estimator can be done as a result of prop score estimates, so can fit back into the GRF framekwork.
firststage_ATE
weighted.mean(firststage.pred$predictions[complier.score > 0],
    1 / complier.score[complier.score > 0])






#! See Estimator 2 Test-and-Select of Hazard (2022), which proposes a
#! cross-fitting approach to this.
#! Think about using the forest to find partitions G (defined by Hazard)
#! The valid approach runs a hypothesis test at the end nodes of trees in forest
# Consider writing up the above first, and then asking Yiqi for her opinion.
#! In particular Athey Wagner (2017 Econometrica) go through the example of
#! efficient policy learning, and classify those who benefit from treatment
#TODO: Implement their approach for those affected by instrument and/or
#TODO: validating IV assumptions.


#! FOUND THE PACKAGE WHICH PREDICTS WHO SHOULD GET THE POLICY:
#! https://github.com/grf-labs/policytree
#! Use this to predict who the compliers are, by the Hazard cross-fitting procedure.

#! Otherwise, I can just do a naive partition of data by X_i (e.g. by naive clustering), from where I can run a local hypothesis test H_0 : \pi (x) = 0 for the region of x decided by the clustering.


################################################################################
## Simple script to show efficiency gain of identifying the compliers:
n <- 10000
Z <- rbinom(n, 1, 0.5)
Q <- rbinom(n, 1, 0.5)
W <- Q * Z
tau <- rep(1, n / 2) %>% c(rep(0, n / 2))
Y <- tau * W + Q + rnorm(n)
# First stage:
lm(W ~ 1 + Z) %>% summary()
# Naive OLS
lm(Y ~ 1 + W) %>% summary()
# IV estimate of LATE
ivreg::ivreg(Y ~ W | Z) %>% summary()
# IV estimate of LATE for compliers -> more efficient, requires oracle knowledge
lm(W[Q == 1] ~ 1 + Z[Q == 1]) %>% summary()
ivreg::ivreg(Y[Q == 1] ~ W[Q == 1] | Z[Q == 1]) %>% summary()
mean(tau)


# Het TE between compliers + never takers
tau <- ifelse(W == 1, 0.75, 0.25)
Y <- tau * W + Q + rnorm(n)
# First stage:
lm(W ~ 1 + Z) %>% summary()
# IV estimate of LATE
ivreg::ivreg(Y ~ W | Z) %>% summary()
# IV estimate of LATE for compliers -> more efficient, requires oracle knowledge
lm(W[Q == 1] ~ 1 + Z[Q == 1]) %>% summary()
ivreg::ivreg(Y[Q == 1] ~ W[Q == 1] | Z[Q == 1]) %>% summary()
