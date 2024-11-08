#!/usr/bin/R
## Senan Hogan-Hennessy, 04 October 2023
## Simulate the COI estimator, whith simulated data following the assumptions.
set.seed(47)
print(Sys.time())

# Functions for data manipulation and visualisation
library(tidyverse)
# Generalised Random Forests, https://grf-labs.github.io/grf/
library(grf)
# The standard, linear, IV estimator package.
library(ivreg)
# Abadie (2003) IV-weighting estimator, An (2016)
# https://github.com/cran/LARF/blob/master/R/larf.R
library(LARF)
# Causal medation package, Imai Keele Yamamoto (2010)
library(mediation)
# Define number of digits in tables and graphs
digits.no <- 3


################################################################################
## Simulate data.

# Sample size
sample.size <- 10^3
# Number of dimensions for the covariates
covariate.dim <- 10
# Generate matrix of covariates, along which treatment intensity varies
min_X <- 0
max_X <- 1 / covariate.dim
X <- runif(sample.size * covariate.dim, min = min_X, max = max_X)
X <- matrix(X, sample.size, covariate.dim)
# Specifiy the treatment groups, etc.
treatment_index_X <- rowSums(X)
idiosyncracy <- 0 + rnorm(sample.size, mean = 0, sd = 0.2)
# some are non-compliers pi(x) = 0, some are compliers pi(X) = 1
complier.threshold <- 1 / 3
alwaystaker.threshold <- 2 / 3
## Model complaince as the result of the result of a latent index (Vycatil 2002)
D_0 <- D_1 <- rep(0, sample.size)
# Never-takers have D_0 = D_1 = 0
D_1[treatment_index_X + idiosyncracy > complier.threshold] <- 1
# Always takers  have D_0 = D_1 = 1
D_0[treatment_index_X + idiosyncracy > alwaystaker.threshold] <- 1
# Label everyone
firststage.designation <- replicate(sample.size, "", "vector")
firststage.designation[D_1 == 0 & D_0 == 0] <- "Never-taker"
firststage.designation[D_1 == 1 & D_0 == 0] <- "Complier"
firststage.designation[D_1 == 1 & D_0 == 1] <- "Always-taker"
table(firststage.designation) / length(firststage.designation)
# Generate the list of a random instrument.
probZ <- 0.5
Z <- rbinom(sample.size, 1, probZ)
# Generate the list of observed binary treatment
D <- (Z * D_1) + ((1 - Z) * D_0)
# Generate the Y Potential Outcomes (POs), works with non-linearities.
outcome_index_X <- 5 * sin(treatment_index_X * pi)
Y_1_1 <- 5 + 2 * outcome_index_X
Y_1_0 <- 10 - 2 * outcome_index_X
Y_0_1 <- 0 + outcome_index_X
Y_0_0 <- 5 - outcome_index_X
#! TEST: POs are treatment gains, along characterics X
Y_1_1 <- (4 * treatment_index_X)^2
Y_1_0 <- (2 * treatment_index_X)^2
Y_0_1 <- (3 * treatment_index_X)^2
Y_0_0 <- (1 * treatment_index_X)^2
#! TEST: Assume exclusion holds
#Y_1_1 <- Y_0_1
#Y_1_0 <- Y_0_0
# Generate the list of observed outcomes
Y <- (Z * D * Y_1_1) +
    (Z * (1 - D) * Y_1_0) +
    ((1 - Z) * D * Y_0_1) +
    ((1 - Z) * (1 - D) * Y_0_0)
# Put these data to a coherent dataframe.
combined.data <- data.frame(
    # Observed data
    Z, D, Y, X,
    # Unobserved, potential outcomes and compliance.
    D_1, D_0,
    Y_1_1, Y_1_0, Y_0_1, Y_0_0,
    treatment_index_X,
    idiosyncracy,
    outcome_index_X,
    firststage.designation) %>%
    tibble()


################################################################################
## Investigate the distributions.

# Plot compliance rate
plot(treatment_index_X, D_1 - D_0)
hist(D_1 - D_0)
# Plot the POs
combined.data %>%
    ggplot(aes(x = treatment_index_X + idiosyncracy)) +
    geom_point(aes(y = Y_1_1, colour = "Y(1, 1)")) +
    geom_point(aes(y = Y_1_0, colour = "Y(1, 0)")) +
    geom_point(aes(y = Y_0_1, colour = "Y(0, 1)")) +
    geom_point(aes(y = Y_0_0, colour = "Y(0, 0)")) +
    theme_bw() +
    scale_x_continuous(name = "Latent Treatment Index, V(X; Z)") +
    scale_y_continuous(name = "") +
    ggtitle("Y(Z, D)") +
    theme(plot.title = element_text(size = rel(1)),
        plot.margin = unit(c(0.5, 0, 0, 0), "mm"))
# Total effect Z -> Y, E[Y | Z = 1] - E[Y | Z = 0] = E[Y(1,D(1)) - Y(0,D(0))]
total_effect <- (Y_1_1 - Y_0_0) * (firststage.designation == "Complier") +
    (Y_1_1 - Y_0_1) * (firststage.designation == "Always-taker") +
    (Y_1_0 - Y_0_0) * (firststage.designation == "Never-taker")
print(mean(total_effect))
# First-stage effect Z -> D, E[D | Z = 1] - E[D | Z = 0] = E[D(1) - D(0)]
firststage_effect <- D_1 - D_0
print(mean(firststage_effect))
# Direct effect D -> Y, E[ Y(Z,1) - Y(Z,0) | Z], for each Z = 0, 1
print(mean(Y_1_1 - Y_1_0))
print(mean(Y_0_1 - Y_0_0))
# Mean AME, E[ E[ Y(Z,1) - Y(Z,0) | Z] ]
mechanism_effect <- Z * (Y_1_1 - Y_1_0) + (1 - Z) * (Y_0_1 - Y_0_0)
print(mean(mechanism_effect))
# Estimand: Mean LAME, local to compliers
lame <- mechanism_effect[firststage.designation == "Complier"]
print(mean(lame))
# Direct effect Z -> Y, E[ Y(1,D) - Y(0,D) | D], for each D = 0, 1
print(mean(Y_1_1 - Y_0_1))
print(mean(Y_1_0 - Y_0_0))
# Estimand: Mean ADE, E[ E[ Y(1,D) - Y(0,D) | D] ]
direct_effect <- D * (Y_1_1 - Y_0_1) + (1 - D) * (Y_1_0 - Y_0_0)
print(mean(direct_effect))
# The Wald estimand, when there is exclusion (Theorem 1., Hogan-Hennessy 2024)
wald_estimate <- (mean(Y[Z == 1]) - mean(Y[Z == 0])) / (
    mean(D[Z == 1]) - mean(D[Z == 0]))
wald_estimand <-
    # Local complier mechanism effect
    mean((Z * (Y_1_1 - Y_1_0) + (1 - Z) * (Y_0_1 - Y_0_0))[
        firststage.designation == "Complier"]) +
    # Local complier direct effect
    mean((Z * (Y_1_0 - Y_0_0) + (1 - Z) * (Y_1_1 - Y_0_1))[
        firststage.designation == "Complier"]) +
    # Local AT direct effect, weighted by Pr(AT) / Pr(complier)
    (mean(firststage.designation == "Always-taker") /
        mean(firststage.designation == "Complier")) * mean(
            (Y_1_1 - Y_0_1)[firststage.designation == "Always-taker"]) +
    # Local NT direct effect, weighted by Pr(NT) / Pr(complier)
    (mean(firststage.designation == "Never-taker") /
        mean(firststage.designation == "Complier")) * mean(
            (Y_1_0 - Y_0_0)[firststage.designation == "Never-taker"])


################################################################################
## Standard IV estimation.

# Linear OLS between educ + income
lm(D ~ 1 + Z) %>% summary()
# The reduced from (Z -> Y, without controlling for D)
lm(Y ~ 1 + Z) %>% summary()
print(mean(total_effect))
# Naive OLS
lm(Y ~ 1 + D) %>% summary()
# Naive 2SLS (without concern for exclusion)
ivreg::ivreg(Y ~ 1 + D | 1 + Z) %>% summary()
ivreg::ivreg(Y ~ 1 + D + X | 1 + Z + X) %>% summary()
print(c(wald_estimate, wald_estimand))
# Abadie (2003) IV-weighting estimator
abadie.reg <- larf(Y ~ 1 + X, treatment = D, instrument = Z, method = "LS")
summary(abadie.reg)
# Estimate the weighting estimator by hand, following Kasy (2016) slides
hat_probZ <- glm(Z ~ 1, family = binomial(probit))$fitted
k_1 <- D * ((Z - hat_probZ) / ((1 - hat_probZ) * hat_probZ))
k_0 <- (1 - D) * (((1 - Z) - (1 - hat_probZ)) / ((1 - hat_probZ) * hat_probZ))
k <- k_1 * hat_probZ + k_0 * (1 - hat_probZ)
# SHow the estiamte, vs the Wald estimand.
print(mean(k_1 * Y) / mean(k_1) - mean(k_0 * Y) / mean(k_0))
print(wald_estimand)
# Hand-code Weighted least Squares solution (Abadie 2003, Section 4.2.1)
hat_probZ <- glm(Z ~ 1 + X, family = binomial(probit))$fitted
k_1 <- D * ((Z - hat_probZ) / ((1 - hat_probZ) * hat_probZ))
k_0 <- (1 - D) * (((1 - Z) - (1 - hat_probZ)) / ((1 - hat_probZ) * hat_probZ))
k <- k_1 * hat_probZ + k_0 * (1 - hat_probZ)
solve(t(cbind(D, 1, X) * k) %*% cbind(D, 1, X)) %*% t(cbind(D, 1, X) * k) %*% as.matrix(Y)
print(wald_estimand)


################################################################################
## Standard IV estimation.

# Linear OLS between educ + income
lm(D ~ 1 + Z) %>% summary()
# The reduced from (Z -> Y, without controlling for D)
lm(Y ~ 1 + Z) %>% summary()
print(mean(total_effect))
# Naive OLS
lm(Y ~ 1 + D) %>% summary()
# Naive 2SLS (without concern for exclusion)
ivreg::ivreg(Y ~ 1 + D | 1 + Z) %>% summary()
ivreg::ivreg(Y ~ 1 + D + X | 1 + Z + X) %>% summary()
print(c(wald_estimate, wald_estimand))
# Abadie (2003) IV-weighting estimator
abadie.reg <- larf(Y ~ 1 + X, treatment = D, instrument = Z, method = "LS")
summary(abadie.reg)
# Estimate the weighting estimator by hand, following Kasy (2016) slides
hat_probZ <- glm(Z ~ 1, family = binomial(probit))$fitted
k_1 <- D * ((Z - hat_probZ) / ((1 - hat_probZ) * hat_probZ))
k_0 <- (1 - D) * (((1 - Z) - (1 - hat_probZ)) / ((1 - hat_probZ) * hat_probZ))
k <- k_1 * hat_probZ + k_0 * (1 - hat_probZ)
# SHow the estiamte, vs the Wald estimand.
print(mean(k_1 * Y) / mean(k_1) - mean(k_0 * Y) / mean(k_0))
print(wald_estimand)
# Hand-code Weighted least Squares solution (Abadie 2003 Section 4.2.1)
hat_probZ <- glm(Z ~ 1 + X, family = binomial(probit))$fitted
k_1 <- D * ((Z - hat_probZ) / ((1 - hat_probZ) * hat_probZ))
k_0 <- (1 - D) * (((1 - Z) - (1 - hat_probZ)) / ((1 - hat_probZ) * hat_probZ))
k <- k_1 * hat_probZ + k_0 * (1 - hat_probZ)
solve(t(cbind(D, 1, X) * k) %*% cbind(D, 1, X)) %*% t(cbind(D, 1, X) * k
    ) %*% as.matrix(Y)
print(wald_estimand)


################################################################################
## New CCI weighting approach.

# Import the pre-defined CCI estimator function, cci.est()
source("cci-define.R")

# Apply the estimator on simulated data.
count.trees <- 2000
cci.weights <- cci_weights.est(Y, D, Z, X, count.trees = count.trees)
cci.est <- cci_point.est(Y, D, Z, X, cci.weights)

# Show how close the estimates are, compared to simulated estimands.
print(c("Total Effect:", cci.est["total_effect"], mean(total_effect)))
print(c("LAME:", cci.est["LAME"],                 mean(lame)))
print(c("AME:", cci.est["mechanism_effect"],      mean(mechanism_effect)))
print(c("ADE:", cci.est["direct_effect"],         mean(direct_effect)))
print(c("AIE:", cci.est["indirect_effect"],       mean(total_effect - direct_effect)))
# Show the partially identified bounds
print(c("AME bounds:",
    cci.est["mechanism_lower"], cci.est["mechanism_effect"],
    "Real value:", mean(mechanism_effect)))
# SHow the naive OLS + IV terms
print(c("regular OLS:", cci.est["ols_est"]))
print(c("regular IV:", cci.est["iv_est"]))

# Do the bounds cover when X under-specified?
cci_bounds.weights <- cci_weights.est(Y, D, Z, X[, 1:2], count.trees = count.trees)
cci_bounds.est <- cci_point.est(Y, D, Z, X[, 1:2], cci.weights)
print(c("AME bounds:",
    cci_bounds.est["mechanism_lower"], cci_bounds.est["mechanism_effect"],
    "Real value:", mean(mechanism_effect)))

# Calculate the boostrap confidence intervals
cci_bootstrap.data <- cci_bootstrap.est(Y, D, Z, X, count.trees = count.trees,
    boot.count = 10)

# Plot the estimates, separately
library(ggridges)
# Get the mean of boot estimates.
mean_estimate.data <- cci_bootstrap.data %>%
    dplyr::select(ols_est, iv_est, LAME) %>%
    pivot_longer(everything(),
        names_to = "parameter", values_to = "estimate") %>%
    group_by(parameter) %>%
    summarise(mean_estimate = mean(estimate)) %>%
    ungroup()
# Plot the distribution of bootestimates
cci_bootstrap.data %>%
    dplyr::select(ols_est, iv_est, LAME) %>%
    pivot_longer(everything(),
        names_to = "parameter", values_to = "estimate") %>%
    ggplot(aes(x = estimate, y = parameter,
        fill = parameter, colour = parameter)) +
    geom_density_ridges2(alpha = 0.5) +
    geom_vline(data = mean_estimate.data,
        aes(xintercept = mean_estimate, colour = parameter),
        linetype = "dashed") +
    theme_bw() +
    scale_x_continuous(name = "Estimate") +
    scale_y_discrete(name = "",
        limits = c("LAME", "iv_est", "ols_est"),
        breaks = c("ols_est", "iv_est", "LAME"),
        labels = c("OLS", "Naive IV", "CCI")) +
    theme(
        #plot.title = element_text(size = rel(1)),
        plot.margin = unit(c(0.5, 0, 0, 0), "mm"),
        legend.position = "none",
        legend.margin = margin(t = -10))
# Show the 95\% percentiles of the relevant parameters
cci_bootstrap.data %>% pull(ols_est) %>% quantile(probs = c(0.025, 0.975)) %>% print()
cci_bootstrap.data %>% pull(iv_est) %>% quantile(probs = c(0.025, 0.975)) %>% print()
cci_bootstrap.data %>% pull(LAME) %>% quantile(probs = c(0.025, 0.975)) %>% print()


#! Untested from here:
quit("no")

#! Test, estimate AME by separate causal forests.
# Estimate Pr(Z = 1 | X) instrument propensity
instpropensity.forest <- probability_forest(X, as.factor(Z),
    num.trees = count.trees)
# Estimate Pr(complier | X) = E[D | X, Z = 1] - E[D | X, Z = 0] by honest GRF
complier.forest <- causal_forest(X, D, Z, num.trees = count.trees)
# Estimate E[Y | X, Z = D = 1] - E[Y | X, Z = 1, D = 0] by causal forest
mechanism_Z1.forest <- causal_forest(X[Z == 1, ], Y[Z == 1], D[Z == 1],
    #Y.hat = Y[Z == 1], W.hat = D[Z == 1],
    num.trees = count.trees)
# Estimate E[Y | X, Z = 0, D = 1] - E[Y | X, Z = D = 0] by causal forest
mechanism_Z0.forest <- causal_forest(X[Z == 0, ], Y[Z == 0], D[Z == 0],
    #Y.hat = Y[Z == 0], W.hat = D[Z == 0],
    num.trees = count.trees)
# Predict across X
inst.propensity <- predict(instpropensity.forest, X)$predictions[, 2]
complier.score <- predict(complier.forest, X)$predictions
mechanism_Z1.prediction <- predict(mechanism_Z1.forest, newdata = X)$predictions
mechanism_Z0.prediction <- predict(mechanism_Z0.forest, newdata = X)$predictions
# Get the treatment effect
mechanism_effect.est <- inst.propensity * mechanism_Z1.prediction +
    (1 - inst.propensity) * mechanism_Z0.prediction


average_treatment_effect(mechanism_Z1.forest)
mean(Y_1_1 - Y_1_0)
average_treatment_effect(mechanism_Z0.forest)
mean(Y_0_1 - Y_0_0)


print(c("AME:", mean(mechanism_effect.est), mean(mechanism_effect)))
print(c("LAME:", weighted.mean(mechanism_effect.est, complier.score), mean(lame)))


plot(lame, mechanism_effect.est[firststage.designation == "Complier"])


data.frame(treatment_index_X, complier.score, firststage.designation) %>%
    ggplot(aes(x = treatment_index_X, y = complier.score,
        colour = firststage.designation)) +
    geom_point()






################################################################################
## Compare estimates using causal medation package (Imai Keele Yamamoto 2010)

# Define the first stage / medation model
first_stage.reg <- lm(D ~ 1 + Z, data = combined.data)
# Define the second stage
second_stage.reg <- lm(Y ~ 1 + Z + D, data = combined.data)
# Estimate the mechanism model.
mechanism.reg <- mediate(first_stage.reg, second_stage.reg,
    treat = "Z", mediator = "D", robustSE = FALSE, sims = 100)
# Show the mechanism estimates -> they are wrong when D is not independent.
summary(mechanism.reg)
print(c("Total Effect:", mean(total_effect)))
print(c("ADE:", mean(direct_effect)))
print(c("AIE:", mean(total_effect - direct_effect)))


# lm(medv ~ poly(lstat, 2, raw = TRUE), data = train.data)
# lm(Y ~ 1 + poly(X, 3, raw = TRUE), subset = (Z == 0 & D == 0))

################################################################################
##! Test: Hand adjustments to the functions


##! TEST: Linear specifivation.
A <- 1
pi <- 1
beta <- 2
g <- 0
Z <- rbinom(sample.size, 1, 0.25)
Q <- rbinom(sample.size, 1, 0.25)
D <- as.integer(Q + Z)
Y <- A + beta * D + Q + g * Z + rnorm(sample.size, mean = 0, sd = 0.1)
lm(D ~ 1 + Z)
lm(Y ~ 1 + D)
ivreg::ivreg(Y ~ 1 + D | 1 + Z) %>% summary()

# Apply the estimator on simulated data.
cci.weights <- cci_weights.est(Y, D, Z, X, num.trees = 500)
cci.est <- cci_direct.est(Y, D, Z, X, cci.weights)

# SHow how close the estimates are, compared to simulated estimands.
print(c("Total Effect:", cci.est["total_effect"]))
lm(Y ~ 1 + Z) %>% summary()
print(c("AME:", cci.est["mechanism_effect"]))
print(c("ADE:", cci.est["direct_effect"]))

# I.e., generate data with D = A + pi Z + e, Y = C + B D + g Z + e.  Does it accurately get B, D?

#! Finding: need two-sided non-compliance, which is not satisified in a linear model
#! Can adjust only one side of the AME, under one-sdied non-compliance.
# TODO adjust the linear model to have two-sided non-compliance.
# A sufficient condition for CCI is that the treatment latent index is not a function of anything else (with non-zero mean) except Z, X.  This might be testable(?).


################################################################################
## Computationally test the COI identification.

#! Note: COI holds true if latent_index_X !corr outcome_index_X

# Test whether COI holds.
combined.data %>%
    mutate(X_bin = ntile(treatment_index_X, n = 100)) %>%
    group_by(X_bin, firststage.designation) %>%
    summarise(
        X     = mean(treatment_index_X),
        Z     = mean(Z),
        Y_1_1 = mean(Y_1_1),
        Y_1_0 = mean(Y_1_0),
        Y_0_1 = mean(Y_0_1),
        Y_0_0 = mean(Y_0_0),
        Y     = mean(Y)) %>%
    ggplot(aes(x = X, y = Y_0_0, colour = firststage.designation)) +
    geom_point() +
    geom_line() +
    geom_smooth()

#! Test whether E[(1 - Z)DY | X] = E[Y(0,1) | X, D(1)=D(0)=1] Pr(Z=0, D=1 | X)
test.data <- combined.data %>%
    filter(firststage.designation == "Always-taker") %>%
    mutate(X_bin = ntile(treatment_index_X, n = 100)) %>%
    group_by(X_bin) %>%
    summarise(Y_0_1 = mean(Y_0_1)) %>%
    ungroup()

combined.data %>%
    mutate(X_bin = ntile(treatment_index_X, n = 100)) %>%
    group_by(X_bin) %>%
    summarise(
        outcome = mean((1 - Z) * D * Y),
        ZD     = mean((1 - Z) * D)) %>%
    ungroup() %>%
    left_join(test.data) %>%
    ggplot(aes(x = outcome, y = Y_0_1 * ZD)) +
    geom_point() +
    geom_line() +
    geom_smooth() +
    geom_abline(slope = 1, intercept = 0)

combined.data %>%
    mutate(X_bin = ntile(treatment_index_X, n = 100)) %>%
    group_by(X_bin) %>%
    summarise(
        outcome = mean((1 - Z) * D * Y),
        ZD     = mean((1 - Z) * D)) %>%
    ungroup() %>%
    left_join(test.data) %>%
    transmute(X_bin = X_bin,
        diff = outcome - (Y_0_1 * ZD)) %>% pull(diff) %>% hist()
