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
sample.size <- 10^4
# Number of dimensions for the covariates
covariate.dim <- 10
# Generate matrix of covariates, along which treatment intensity varies
min_X <- 0
max_X <- 1 / covariate.dim
X <- runif(sample.size * covariate.dim, min = min_X, max = max_X)
X <- matrix(X, sample.size, covariate.dim)
# Specifiy the treatment groups, based on covariates + potential outcomes
treatment_index_X <- rowSums(X)
#! TEST: POs are treatment gains, along characterics X
Y_0_0 <- 1 + 10 * treatment_index_X + sin(15 * treatment_index_X)
Y_0_1 <- 15 * treatment_index_X + 2 * cos(15 * treatment_index_X)
Y_1_0 <- 2 + 20 * treatment_index_X + sin(10 * treatment_index_X)
Y_1_1 <- 2.5 + 20 * treatment_index_X + 3 * cos(20 * treatment_index_X)
#! TEST: Assume exclusion holds
#Y_0_1 <- Y_1_1
#Y_0_0 <- Y_1_0
# Model complaince as imperfect selection (i.e., latent index Vycatil 2002)
D_0 <- D_1 <- rep(0, sample.size)
D_1 <- as.integer(Y_1_1 + rnorm(sample.size, 0, 3) >= Y_1_0)
D_0 <- as.integer(Y_0_1 + rnorm(sample.size, 0, 3) >= Y_0_0)
# Replace treatment gains where there would be defiers
#D_0[D_1 == 0 & D_0 == 1] <- 0
# Label compliance
firststage.designation <- replicate(sample.size, "", "vector")
firststage.designation[D_1 == 0 & D_0 == 0] <- "Never-taker"
firststage.designation[D_1 == 1 & D_0 == 1] <- "Always-taker"
firststage.designation[D_1 == 1 & D_0 == 0] <- "Complier"
firststage.designation[D_1 == 0 & D_0 == 1] <- "Defier"
table(firststage.designation) / length(firststage.designation)
# Generate the list of a random instrument.
probZ <- 0.5
Z <- rbinom(sample.size, 1, probZ)
# Generate the list of observed binary treatment
D <- (Z * D_1) + ((1 - Z) * D_0)
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
    firststage.designation) %>%
    tibble()


################################################################################
## Investigate the distributions.

# Plot compliance rate
plot(treatment_index_X, D_1 - D_0)
hist(D_1 - D_0)
# Plot the POs
combined.data %>%
    ggplot(aes(x = treatment_index_X)) +
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

print(wald_estimate)
summary(lm(Y ~ 1 + D + X, weights = k))

# Hand-code Weighted least Squares solution (Abadie 2003 Section 4.2.1)
hat_probZ <- glm(Z ~ 1 + poly(X, 3), family = binomial(probit))$fitted
k_1 <- D * ((Z - hat_probZ) / ((1 - hat_probZ) * hat_probZ))
k_0 <- (1 - D) * (((1 - Z) - (1 - hat_probZ)) / ((1 - hat_probZ) * hat_probZ))
k <- k_1 * hat_probZ + k_0 * (1 - hat_probZ)
solve(t(cbind(D, 1, X) * k) %*% cbind(D, 1, X)) %*% t(cbind(D, 1, X) * k) %*% as.matrix(Y)
print(wald_estimate)
print(wald_estimand)

#! Test: is E[ Y(0, 0) | D(1) = 1, D(0) = 0] identified?
mean(Y_0_0[firststage.designation == "Complier"])
mean(Y[D == 0 & firststage.designation == "Complier"])
mean((Y * (1 - D) / (1 - hat_probZ))[firststage.designation == "Complier"])
mean(Y * k * (1 - D) / (1 - hat_probZ)) / mean(k)
mean(k_0 * Y) / mean(k_0)

#! Test: is E[ Y(1, 1) | D(1) = 1, D(0) = 0] identified?
mean(Y_1_1[firststage.designation == "Complier"])
mean(Y[D == 1 & firststage.designation == "Complier"])
mean((Y * D / hat_probZ)[firststage.designation == "Complier"])
mean((k * Y * D / hat_probZ)) / mean(k)
mean(k_1 * Y) / mean(k_1)


#!RESULT: E[ k * g(Y) ] / mean(k) != E[ g(Y) | complier]
#!RESULT: the proof breaks down as
#!RESULT: E[ g(Y) | D(1) = D(0) = 1] x Pr(D(1) = D(0) = 1) is not equal to
#!RESULT: E[ g(Y) | D = 1, Z = 0] x Pr(D = 1 | Z = 0)
# See this for the standard case,
mean(k * Y) / mean(k)
mean(Y[firststage.designation == "Complier"])
# This would have been true for g(X), but Y has manipulation along the effect of Z
#!Thoughts: not sure why this is, and have little time to figure it out.

mean(k * Y) / mean(k)
mean(Y[firststage.designation == "Complier"])
mean((Z * Y_1_1 + (1 - Z) * Y_0_0)[firststage.designation == "Complier"])


#! Test, what do the k weights estimate now?
mean(D * (1 - Z) * Y) / mean(D * (1 - Z))
mean(Y[firststage.designation == "Always-taker"])

mean(k * Y) / mean(k)
mean(Y[firststage.designation == "Complier"])
mean(k_0 * Y) / mean(k_0)
mean(Y_0_0[firststage.designation == "Complier"])
mean(k_1 * Y) / mean(k_1)
mean(Y_1_1[firststage.designation == "Complier"])
mean(k_1 * Y) / mean(k_1) - mean(k_0 * Y) / mean(k_0)
mean((Y_1_1 - Y_0_0)[firststage.designation == "Complier"])


# WHat happens in Abadie (2002) estimator?
(mean((Y * D)[Z == 1]) - mean((Y * D)[Z == 0])) / (
    mean(D[Z == 1]) - mean(D[Z == 0]))
mean(Y_1_1[firststage.designation == "Complier"])

(mean((Y * (1 - D))[Z == 1]) - mean((Y * (1 - D))[Z == 0])) / (
    mean((1 - D)[Z == 1]) - mean((1 - D)[Z == 0]))
mean(Y_0_0[firststage.designation == "Complier"])

