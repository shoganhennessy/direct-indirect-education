# Show the date:
print(format(Sys.time(), "%H:%M %Z %A, %d %B %Y"))
# Set up the R environment
set.seed(47)
# Load tidyverse, for dealing with data and plotting.
suppressPackageStartupMessages(
    library("tidyverse", quietly = TRUE))
# Load the library for linear causal mediation (Imai, Keele, Yamamoto 2010).
suppressPackageStartupMessages(
    library("mediation", quietly = TRUE))
# Package to annotate ggplot with TeX maths
suppressPackageStartupMessages(
    library("latex2exp", quietly = TRUE))
# Library for better colour choice.
suppressPackageStartupMessages(
    library("ggthemes", quietly = TRUE))
# Set the options for the plot sizes
fig.height <- 4.5
fig.width <- 1.25 * fig.height
options(repr.plot.width = fig.width, repr.plot.height = fig.height,
    repr.plot.res = 175)


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
# Specification: Generate outcomes as a sum of observed, mu(X), + unobserved, U
treatment_index_X <- rowSums(X)
mu0_Z0_X <- rep(10, sample.size)
mu1_Z0_X <- 0.9 * mu0_Z0_X + 25 * as.integer(treatment_index_X >= 0.55)
mu0_Z1_X <- mu0_Z0_X * 1.2
mu1_Z1_X <- 0.9 * mu0_Z1_X + 50 * as.integer(treatment_index_X >= 0.45)
# To aid theory on selection-into-educ, suppose educ gains >= 0
#mu1_Z0_X[mu1_Z0_X < mu0_Z0_X] <- mu0_Z0_X
#mu1_Z1_X[mu1_Z1_X < mu0_Z1_X] <- mu0_Z1_X
# Unobserved selection, between D = 1 and D = 0
# Assume: U_1, U_0 ~ BivarNormal(rho, mu_0 = mu_1 = 0, sigma_0, sigma_1)
rho <- 0.75
sigma_0 <- 1
sigma_1 <- 2 * sigma_0
U_both <- MASS::mvrnorm(
    n = sample.size,
    mu = c(0, 0),
    Sigma = matrix(c(
        sigma_0^2, rho * sigma_0 * sigma_1,
        rho * sigma_0 * sigma_1, sigma_1^2
    ), ncol = 2)
)
U_0 <- U_both[, 1]
U_1 <- U_both[, 2]
print(paste0("sigma_0 = ", sigma_0, ", sigma_1 = ", sigma_1))
# Compose potential outcomes from the observed + unobserved factors.
# Y_i(Z, D) = mu_D(Z; X_i) + U_D
Y_0_0 <- mu0_Z0_X + U_0
Y_0_1 <- mu1_Z0_X + U_1
Y_1_0 <- mu0_Z1_X + U_0
Y_1_1 <- mu1_Z1_X + U_1
# Model compliance as pure selection --- Roy gains as latent index/prop score.
# D(Z) = 1{ Y_i(Z, D = 0)    <= Y(Z, D_i = 1) }
#      = 1{ mu_0(Z; X_i) + U_0 <= mu_1(Z; X_i) + U_1}
D_0 <- D_1 <- replicate(0, sample.size)
D_1 <- as.integer(Y_1_0 <= Y_1_1)
D_0 <- as.integer(Y_0_0 <= Y_0_1)
# Label compliance
firststage.designation <- replicate(sample.size, "", "vector")
firststage.designation[D_1 == 0 & D_0 == 0] <- "Never-taker"
firststage.designation[D_1 == 1 & D_0 == 1] <- "Always-taker"
firststage.designation[D_1 == 1 & D_0 == 0] <- "Complier"
firststage.designation[D_1 == 0 & D_0 == 1] <- "Defier"
print(table(firststage.designation) / length(firststage.designation))

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

# Plot the POs
combined.data %>%
    #filter(firststage.designation == "Complier") %>%
    ggplot(aes(x = treatment_index_X)) +
    geom_point(aes(y = Y_1_1, colour = "Y(1, 1)")) +
    geom_point(aes(y = Y_1_0, colour = "Y(1, 0)")) +
    geom_point(aes(y = Y_0_1, colour = "Y(0, 1)")) +
    geom_point(aes(y = Y_0_0, colour = "Y(0, 0)")) +
    theme_bw() +
    scale_x_continuous(name = "Latent Treatment Index, V(X)") +
    scale_y_continuous(name = "") +
    ggtitle("Y(Z, D)") +
    theme(
        plot.title = element_text(size = rel(1)),
        plot.margin = unit(c(0.5, 0, 0, 0), "mm"))

#! Testing to get 50-50 indirect + direct effects.
# Get the theoretical total effect/reduced form/ITT
total_effect <- (Y_1_1 - Y_0_0) * (firststage.designation == "Complier") +
    (Y_1_1 - Y_0_1) * (firststage.designation == "Always-taker") +
    (Y_1_0 - Y_0_0) * (firststage.designation == "Never-taker")
average_total_effect <- mean(total_effect)
print(paste("Average total effect:", round(average_total_effect, digits.no)))
# Get the theoretical indirect effect
indirect_effect <- (Z * (Y_1_1 - Y_1_0) + (1 - Z) * (Y_0_1 - Y_0_0)) * (
        firststage.designation == "Complier") +
    (0) * (firststage.designation == "Always-taker") +
    (0) * (firststage.designation == "Never-taker")
average_indirect_effect <- mean(indirect_effect)
print(paste("Average indirect effect:",
    round(average_indirect_effect, digits.no)))

# Plot the compliance
combined.data %>%
    ggplot(aes(x = treatment_index_X)) +
    # geom_point(aes(y = D_1, colour = "D(Z = 1)")) +
    # geom_point(aes(y = D_0, colour = "D(Z = 0)")) +
    geom_point(aes(y = D_1 - D_0, colour = "red")) +
    geom_smooth(aes(y = D_1 - D_0, colour = "red"), se = FALSE) +
    theme_bw() +
    scale_x_continuous(name = "Latent Treatment Index, V(X)") +
    scale_y_continuous(name = "") +
    ggtitle("D(1) - D(0), compliance across X_i") +
    theme(
        plot.title = element_text(size = rel(1)),
        plot.margin = unit(c(0.5, 0, 0, 0), "mm"))


################################################################################
## First-stage compliance variance decomposition.

# True value of total variance, Var(D_1 - D_0), and simple estimator
var(D_1 - D_0)
p <- mean(D[Z == 1]) - mean(D[Z == 0])
print(p * (1 - p))
# Simple bootstrap appropriate for the SEs, thanks to the delta method (to-do)

# Import the Sanchez-Bacon (2013) package to estimate Var( E[D(1) - D(0) | X])
source(file.path("r-vcate", "base", "estimation_fns.R"))
source(file.path("r-vcate", "base", "wrapper_fns.R"))
# Estimate the first-stage VCATE, conditional on X.
firststage.vcate <- estimate_vcate(
    y = D,
    d = Z,
    x_mat = X,
    px = rep(0.5, nrow(X)), # Pr(Z = 1 | X) = 0.5, for all X
    clustervar = seq_len(nrow(X)), # No clusters
    num_folds = 4,
    family.choice = "binomial", # Logit specification, outcome D binary
    interpret_effect_size = FALSE,
    num_splits = 4,
    alpha_target = 0.1,
    acc = 0.001,
    seed = 47,
    messages = TRUE
)
# Show the estimates for first-stage VCATE.
print(firststage.vcate$ci)
print(firststage.vcate$median_vtauhat)

## Compare to the true value of population compliance explained by X_i
# 1. Notice that we can rewrite D_i(1) as
# D_1 = 1{mu0_Z1_X + U_0 <= mu1_Z1_X + U_1} = 1{U_0-U_1 <= mu1_Z1_X - mu0_Z1_X}
# Since (U_0 - U_1) is normally distributed, then the conditional mean given X
prob_D1_X <- pnorm((mu1_Z1_X - mu0_Z1_X) / sd(U_0 - U_1))
# 2. By similar logic for D_i(0)
# D_0 = 1{mu0_Z0_X + U_0 <= mu1_Z0_X + U_1} = 1{U_0-U_1 <= mu1_Z0_X - mu0_Z0_X}
prob_D0_X <- pnorm((mu1_Z0_X - mu0_Z0_X) / sd(U_0 - U_1))
cate_d_x <- prob_D1_X - prob_D0_X
print(var(cate_d_x))
print(firststage.vcate$ci)
print(firststage.vcate$median_vtauhat)

## Compare to the true value of population mediation. explained by X_i
D_X <- Z * (mu1_Z1_X - mu0_Z1_X) + (1 - Z) * (mu1_Z0_X - mu0_Z0_X)
prob_D_X <- pnorm(D_X / sd(U_0 - U_1))
print(var(prob_D_X))

# Calculate the estimated R_U^2 value
R2_U.est <- 1 - firststage.vcate$median_vtauhat / (p * (1 - p))
print(R2_U.est)
print(1 - var(cate_d_x) / var(D_1 - D_0))


################################################################################
## First-stage compliance decomposition.

# True value of variance of D
print(var(D))
D_X <- Z * (mu1_Z1_X - mu0_Z1_X) + (1 - Z) * (mu1_Z0_X - mu0_Z0_X)
prob_D_X <- pnorm(D_X / sd(U_0 - U_1))
var_prob_D_X <- var(prob_D_X)

# Estimate E[D | Z, X]
firststage.forest <- boosted_regression_forest(cbind(Z, X), D)
firststage.est <- predict(firststage.forest)$predictions
# Show point estimate of Var(E[D | Z, X]), without accounting prediction error.
print(var_prob_D_X)
print(var(firststage.est))
R2_U.est <- 1 - var(firststage.est) / var(D)
print(R2_U.est)


################################################################################
## Linear model of Causal medation package (Imai Keele Yamamoto 2010)

# Define the first stage / medation model
first_stage.reg <- lm(D ~ 1 + Z + X, data = combined.data)
# Total effect
total_stage.reg <- lm(Y ~ 1 + Z + X, data = combined.data)
# Define the second stage
second_stage.reg <- lm(Y ~ 1 + Z + D + Z:D + X, data = combined.data)
# Estimate the mechanism model.
mediation.reg <- mediate(first_stage.reg, second_stage.reg,
    treat = "Z", mediator = "D",
    robustSE = FALSE, sims = 5)
summary(mediation.reg)

# Non-parametric estimation:
# DML estimates of the direct and indirect effects.
library(causalweight)
mediation.dml <- medDML(y = Y, d = Z, m = D, x = X)
print(round(mediation.dml$results, digits.no))

## Show the system estimates.
# First stage Z -> D(Z), E[D | Z = 1] - E[D | Z = 0] = E[D(1) - D(0)]
print(mean(D_1 - D_0))
print(coef(first_stage.reg)["Z"])
# Get the theoretical total effect/reduced form/ITT
total_effect <- (Y_1_1 - Y_0_0) * (firststage.designation == "Complier") +
    (Y_1_1 - Y_0_1) * (firststage.designation == "Always-taker") +
    (Y_1_0 - Y_0_0) * (firststage.designation == "Never-taker")
average_total_effect <- mean(total_effect)
print(paste("Average total effect:", round(average_total_effect, digits.no)))
print(coef(total_stage.reg)["Z"])
# Get the theoretical indirect effect
indirect_effect <- (Z * (Y_1_1 - Y_1_0) + (1 - Z) * (Y_0_1 - Y_0_0)) * (
        firststage.designation == "Complier") +
    (0) * (firststage.designation == "Always-taker") +
    (0) * (firststage.designation == "Never-taker")
average_indirect_effect <- mean(indirect_effect)
print(paste("Average indirect effect:",
    round(average_indirect_effect, digits.no)))
print(mediation.reg$d.avg)
print((coef(second_stage.reg)["D"] + 
    coef(first_stage.reg)["Z"] * coef(second_stage.reg)["Z:D"]) * (
        coef(first_stage.reg)["Z"]))
# Get the theoretical direct effect
direct_effect <- (Z * (Y_1_1 - Y_0_1) + (1 - Z) * (Y_1_0 - Y_0_0)) * (
        firststage.designation == "Complier") +
    (Y_1_1 - Y_0_1) * (firststage.designation == "Always-taker") +
    (Y_1_0 - Y_0_0) * (firststage.designation == "Never-taker")
average_direct_effect <- mean(direct_effect)
print(paste("Average direct effect:", round(average_direct_effect, digits.no)))
print(mediation.reg$z.avg)

# Show all estimates together.
summary(mediation.reg)



################################################################################
## Regular matching estimators (observing X) get it correct.

# Non-parametric matching estimator.
library(MatchIt)
# https://cran.r-project.org/web/packages/MatchIt/vignettes/matching-methods.html
# https://cran.r-project.org/web/packages/MatchIt/vignettes/estimating-effects.html
matching_firststage.est <- matchit(
    D ~ 1 + Z +
        X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10,
    method = "nearest",
    distance = "mahalanobis",
    data = combined.data
)
matching.data <- match.data(matching_firststage.est)
matching_secondstage.est <- lm(Y ~ 1 + D + Z,
    weights = weights,
    data = matching.data
)
marginaleffects::avg_comparisons(matching_secondstage.est,
    variables = "D",
    vcov = ~subclass,
    newdata = matching.data,
    wts = "weights"
)
print(mean((Z[D == 1] * (Y_1_1[D == 1] - Y_1_0[D == 1]) + (1 - Z[D == 1]) * (Y_0_1[D == 1] - Y_0_0[D == 1]))))
print(mean(indirect_effect[D == 1]) / mean((D_1 - D_0)[D == 1]))
# Unknown which TE this is in general
