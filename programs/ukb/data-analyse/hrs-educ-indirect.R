#!/usr/bin/R
## Senan Hogan-Hennessy, 25 September 2024
## Invalid IV of educ gene score for higher ed on earnings.
set.seed(47)
print(Sys.time())

# Functions for data manipulation and visualisation
library(tidyverse)
# Generalised Random Forests, https://grf-labs.github.io/grf/
library(grf)
# The standard, linear, IV estimator package.
library(ivreg)
# Causal medation package, Imai Keele Yamamoto (2010)
library(mediation)
# My custom flavour of Stargazer TeX tables:
# devtools::install_github("shoganhennessy/stargazer")
library(stargazer)
# Define number of digits in tables and graphs
digits.no <- 3
# Size for figures
fig.width <- 15
fig.height <- fig.width * 1 / 2
# List of default colours.
colour.list <- c(
    "#1f77b4", # Blue
    "#d62728", # Red
    "#2ca02c") # Green
# Define folder path for the general data.
data.folder <- file.path("..", "..", "data")
# Graphics output folders
presentation.folder <- file.path("..", "..", "text", "presentation-files")
figures.folder <- file.path("..", "..", "text", "sections", "figures")
tables.folder <- file.path("..", "..", "text", "sections", "tables")


################################################################################
## Load collapsed version of HRS panel data

# Load the pre-cleaned + collapsed HRS panel data.
hrs_panel.data <- data.folder %>%
    file.path("hrs-public", "hrs-panel-collapsed.csv") %>%
    read_csv()


################################################################################
## Set up the variables to study.

# Show variable names.
print(names(hrs_panel.data))

# Filter to non-missings
hrs_panel.data <- hrs_panel.data %>%
    # Filter to non-missing in the input + outcome variables.
    filter(!is.na(genescore_educ_euro),
        !is.na(indiv_edyears),
        !is.na(degree_reported),
        !is.na(indiv_income_real))

# Filter to non-missings
hrs_panel.data <- hrs_panel.data %>%
    mutate(parent_edyears =
        pmax(mother_edyears, father_edyears, na.rm = TRUE)) %>%
    filter(!is.na(parent_edyears)) %>%
    # Put zeros in the missings
    mutate(
        mother_edyears = ifelse(is.na(mother_edyears), 0, mother_edyears),
        father_edyears = ifelse(is.na(father_edyears), 0, mother_edyears))

    
# Select all relevant variables.
train_hrs_panel.data <- hrs_panel.data %>%
    # Select, and define, the relevant variables.
    dplyr::select(
        # Inst Z_i, treatment D_i, outcome Y_i
        genescore_educ_euro, degree_reported, indiv_edyears, indiv_income_real,
        # Demographic controls X_i
        survey_year,
        observed_count,
        hhres,
        #child,
        gender_female,
        indiv_agey,
        parent_edyears,
        mother_edyears,
        father_edyears,
        # Other gene score controls
        genescore_gencog_euro,
        genescore_bmi_euro,
        genescore_waisthips_euro,
        genescore_height_euro,
        genescore_schizo_euro,
        genescore_alcohol_euro,
        genescore_alzh_euro,
        genescore_neurotic_euro,
        genescore_wellbeing_euro,
        genescore_heartcvd_euro,
        genescore_heartmi_euro,
        genescore_bipolar_euro,
        genescore_adhd_euro,
        genescore_crossdisorder_euro,
        genescore_extravert_euro,
        genescore_longevity_euro,
        genescore_antisocial_euro,
        genescore_ocd_euro,
        genescore_mdepressived_euro,
        genescore_ptsd_euro,
        genescore_anxietyfactor_euro,
        genescore_anxietycontrol_euro,
        genescore_agesmoking_euro)

# Filter to only high-school graduates, and college grads
train_hrs_panel.data <- train_hrs_panel.data %>%
    filter(12 <= indiv_edyears & indiv_edyears <= 16)

#! TEST: any missings?
print(c(dim(drop_na(train_hrs_panel.data)), dim(train_hrs_panel.data)))

# Put relevant data into matrix form
Z <- train_hrs_panel.data %>%
    transmute(genescore_educ_euro = as.integer(0 <= genescore_educ_euro)) %>%
    pull(genescore_educ_euro)
D <- train_hrs_panel.data %>%
    transmute(indiv_edyears = as.integer(indiv_edyears >= 14)) %>%
    pull(indiv_edyears)
X <- train_hrs_panel.data %>%
    dplyr::select(- c(degree_reported, indiv_edyears,
        genescore_educ_euro, indiv_income_real),
        starts_with("genescore"), - genescore_educ_euro) %>%
    #! Consider removing the other gene score controls
    dplyr::select(- starts_with("genescore")) %>%
    as.matrix()
Y <- train_hrs_panel.data %>%
    transmute(Y = log(indiv_income_real)) %>%
    pull(Y)


################################################################################
## Analyse direct + indirect effects.

## Test the exclusion restriction by the Kitagawa (2015) test statistic.
# p < \alpha rejects H_0 : randomisation of Z and/or exclusion restriction
source(file.path("..", "simulations", "kitagawa-2015", "kitagawa-test.R"))
xis <- c(0.01, 0.05, 0.1)
xis <- sqrt(xis * (1 - xis))
# print(Z_validity_test(Y, D, Z, c(0, 1), xis))


################################################################################
## Valculate the R^2_U value among the HRS data.

# Simple estimator for Var(D_1 - D_0)
p <- mean(D[Z==1]) - mean(D[Z==0])
print(p * (1 - p))
# Conditional on other stuff
summary(lm(D ~ 1 + Z + X))
p <- as.numeric(coef(lm(D ~ 1 + Z + X))["Z"])
print(p * (1 - p))

# Import the Sanchez-Bacon method to get Var(D(1) - D(0) | X)
source(file.path("..", "simulations", "r-vcate", "base", "estimation_fns.R"))
source(file.path("..", "simulations", "r-vcate", "base", "wrapper_fns.R"))
# Estimate the first-stage VCATE, conditional on X.
firststage.vcate <- estimate_vcate(
    y            = D,
    d            = Z,
    x_mat        = X,
    px           = rep(0.5, nrow(X)),
    clustervar   = seq_len(nrow(X)),
    num_folds    = 2,
    family.choice = "binomial", # Logit specification, outcome D binary
    interpret_effect_size = FALSE,
    num_splits   = 2,
    alpha_target = 0.1,
    acc          = 0.001,
    seed         = 47,
    messages     = TRUE)

# Show the estimates for first-stage VCATE.
print(firststage.vcate$ci)
print(firststage.vcate$median_vtauhat)
print(p * (1 - p))
# Calculate the R^2_U value
R2_U <- 1 - firststage.vcate$median_vtauhat / (p * (1 - p))
print(R2_U)


################################################################################
## First-stage compliance decomposition.

# True value of variance of D
print(var(D))
# Estimate E[D | Z, X]
firststage.forest <- boosted_regression_forest(cbind(Z, X), D, num.trees = 5000)
firststage.est <- predict(firststage.forest)$predictions
# Show point estimate of Var(E[D | Z, X]), without accounting prediction error.
print(var(firststage.est))
R2_U.est <- 1 - var(firststage.est) / var(D)
print(R2_U.est)
#! Poor predictive power in HRS data.


################################################################################
## Estimate by causal mediation, comparing the results.

# Estimate the total effect
print(summary(lm(Y ~ (1 + Z) * X)))

# Estimate the mediation model, first stage is IV first-stage, second is SI
mediation_nocontrols.reg <- mediate(
    lm(D ~ 1 + Z),     # First-stage
    lm(Y ~ 1 + Z * D), # Second-stage
    treat = "Z", mediator = "D",
    robustSE = TRUE, sims = 500)
# Show the mechanism estimates (note: upper bounds when D +ve unobs selection)
print(summary(mediation_nocontrols.reg))

# Estimate the mediation model, first stage is IV first-stage, second is SI
mediation.reg <- mediate(
    lm(D ~ (1 + Z) * X),     # First-stage
    lm(Y ~ (1 + Z * D) * X), # Second-stage
    treat = "Z", mediator = "D",
    robustSE = TRUE, sims = 500)
# Show the mechanism estimates (note: upper bounds when D +ve unobs selection)
print(summary(mediation.reg))

# DML estimates of the direct and indirect effects.
library(causalweight)
mediation.dml <- medDML(y = Y, d = Z, m = D, x = X)
print(round(mediation.dml$results, 3))

# Estimate the mediation with a CF estimator.
## Calculate the direct + indirect effect point estimates by expectations.
library(sampleSelection)
example.data <- tibble(data.frame(Z = Z, D = D, X = X))
# Get the total effect by prediction
totaleffect.reg <- lm(Y ~ (1 + Z) * X, data = example.data)
mean(predict(totaleffect.reg, newdata = mutate(example.data, Z = 1)) -
    predict(totaleffect.reg, newdata = mutate(example.data, Z = 0)))
# Get the first-stage by prediction
firststage.reg <- probit(D ~ (1 + Z) * X, data = example.data)
mean(predict(firststage.reg, newdata = mutate(example.data, Z = 1), type = "response") -
    predict(firststage.reg, newdata = mutate(example.data, Z = 0), type = "response"))
# Get the second-stage by using a parametric control function.
# Unobserved part to add on, + U_0 + I(D * (U_1 - U_0)), see H-H (2024).
# Or estimated with a Heckman (switching) selection approach.
firststage.eq    <- (D ~ (1 + Z) * X)
outcomeY_Z_0.eq  <- (Y ~ (1 + Z) * X)
outcomeY_Z_1.eq  <- (Y ~ (1 + Z) * X)
selection.reg <- selection(
    firststage.eq,
    list(outcomeY_Z_0.eq, outcomeY_Z_1.eq),
    data = example.data,
    method = "ml")
# Direct Effect, E[ Y_i(1, D) - Y_i(0, D)].
mean(predict(selection.reg, newdata = mutate(example.data, Z = 1)) -
    predict(selection.reg, newdata = mutate(example.data, Z = 0)))
# Indirect effect, E[ Y_i(Z, 1) - Y_i(Z, 0) \times (D_i(1) - D_i(0)) ]
mean((
    predict(firststage.reg, newdata = mutate(example.data, Z = 1), type = "response") -
    predict(firststage.reg, newdata = mutate(example.data, Z = 0), type = "response")) * (
        predict(selection.reg, newdata = mutate(example.data, D = 1))[, 2] -
        predict(selection.reg, newdata = mutate(example.data, D = 0))[, 1]))

# TODO: bootstrap the above.
# TODO: Work out how to add in a control function, 2-step approach (+ boot it).

#! TODO: If I want to augment with a DML estimator, then I will have to code it by hand:
#!       https://docs.doubleml.org/stable/intro/intro.html
#!       https://docs.doubleml.org/stable/guide/models.html#irm-model
