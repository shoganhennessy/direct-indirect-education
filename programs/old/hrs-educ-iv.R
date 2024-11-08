#!/usr/bin/R
## Senan Hogan-Hennessy, 25 September 2023
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
    file.path("hrs-panel-collapsed.csv") %>%
    read_csv()


################################################################################
## Set up the variables to study.

names(hrs_panel.data)
# Filter to non-missings
hrs_panel.data <- hrs_panel.data %>%
    # Filter to non-missing in the input + outcome variables.
    filter(!is.na(genescore_educ_euro),
        !is.na(degree_reported), !is.na(indiv_income_real))
    # Restrict to HS vs BA
    #filter(degree_reported %in% 1:5)

# Select all relevant variables.
train_hrs_panel.data <- hrs_panel.data %>%
    # Select, and define, the relevant variables.
    dplyr::select(
        # Inst Z_i, treatment D_i, outcome Y_i
        genescore_educ_euro, degree_reported, indiv_income_real,
        # Demographic controls X_i
        survey_year,    observed_count, hhres, #child,
        gender_female,  indiv_agey,
        parent_edyears, #mother_edyears, father_edyears,
        # Other gene score controls
        genescore_gencog_euro,
        #genescore_bmi_euro,
        #genescore_waisthips_euro,
        genescore_height_euro,
        #genescore_schizo_euro,
        #genescore_alcohol_euro,
        #genescore_alzh_euro,
        #genescore_neurotic_euro,
        #genescore_wellbeing_euro,
        #genescore_heartcvd_euro,
        #genescore_heartmi_euro,
        genescore_bipolar_euro,
        genescore_adhd_euro,
        genescore_crossdisorder_euro,
        genescore_extravert_euro,
        #genescore_longevity_euro,
        #genescore_antisocial_euro,
        #genescore_ocd_euro,
        #genescore_mdepressived_euro,
        #genescore_ptsd_euro,
        #genescore_anxietyfactor_euro,
        #genescore_anxietycontrol_euro,
        genescore_agesmoking_euro)
        # Other outcomes.
        # pmbmi,          pmhght,         pmwght,
        # child,          cendiv,
        # cogtot,         nsscre,         vocab,          imrc,       dlrc

#! TEST: any missings?
dim(drop_na(train_hrs_panel.data)) == dim(train_hrs_panel.data)

# Put relevant data into matrix form
Z <- train_hrs_panel.data %>%
    transmute(genescore_educ_euro = as.integer(genescore_educ_euro >= 0)) %>%
    pull(genescore_educ_euro)
D <- train_hrs_panel.data %>%
    transmute(degree_reported = as.integer(degree_reported >= 5)) %>%
    pull(degree_reported)
X <- train_hrs_panel.data %>%
    dplyr::select(- c(
        genescore_educ_euro, degree_reported, indiv_income_real),
        - starts_with("genescore")) %>%
    as.matrix()
G <- train_hrs_panel.data %>%
    dplyr::select(starts_with("genescore"), - genescore_educ_euro) %>%
    as.matrix()
Y <- train_hrs_panel.data %>%
    transmute(Y = log(indiv_income_real)) %>%
    pull(Y)


################################################################################
## Analyse the reduced form and first stage: edu gene score -> income+education.

# naive OLS between educ + income
wagediff.reg <- lm(Y ~ 1 + D)
wagediff_controls.reg <- lm(Y ~ 1 + D + X)
wagediff_genecontrols.reg <- lm(Y ~ 1 + D + X + G)
summary(wagediff.reg)
summary(wagediff_controls.reg)
summary(wagediff_genecontrols.reg)


# Linear OLS between educ + gene score
reducedform.reg <- lm(Y ~ 1 + Z)
reducedform_controls.reg <- lm(Y ~ 1 + Z + X)
reducedform_genecontrols.reg <- lm(Y ~ 1 + Z + X + G)
summary(reducedform.reg)
summary(reducedform_controls.reg)
summary(reducedform_genecontrols.reg)

# Output to a LaTeX Tables
stargazer(
    reducedform.reg,
    reducedform_controls.reg,
    reducedform_genecontrols.reg,
    dep.var.caption = "Dependent Variables: Log Earnings",
    add.lines = list(
        c("Gene Score Controls?", "No", "No", "Yes")),
    #omit = "gene",
    omit.stat = c("LL", "ser", "aic", "adj.rsq"),
    star.cutoffs = NA,
    header = FALSE, float = FALSE, no.space = TRUE,
    omit.table.layout = "n", notes.append = FALSE,
    type = "text",
    out = file.path(tables.folder, "reducedform-reg.tex"))

# Linear OLS between educ + gene score
firststage.reg <- lm(D ~ 1 + Z)
firststage_controls.reg <- lm(D ~ 1 + Z + X)
firststage_genecontrols.reg <- lm(D ~ 1 + Z + X + G)
summary(firststage.reg)
summary(firststage_controls.reg)
summary(firststage_genecontrols.reg)
# Output to a LaTeX Tables
stargazer(
    firststage.reg,
    firststage_controls.reg,
    firststage_genecontrols.reg,
    dep.var.caption = "Dependent Variables: Attended university",
    add.lines = list(
        c("Gene Score Controls?", "No", "No", "Yes")),
    #omit = "gene",
    omit.stat = c("LL", "ser", "aic", "adj.rsq"),
    star.cutoffs = NA,
    header = FALSE, float = FALSE, no.space = TRUE,
    omit.table.layout = "n", notes.append = FALSE,
    type = "text",
    out = file.path(tables.folder, "firststage-reg.tex"))

# Linear OLS naive IV
secondstage.reg <- ivreg(Y ~ 1 + D | 1 + Z)
secondstage_controls.reg <- ivreg(Y ~ 1 + D + X | 1 + Z + X)
secondstage_genecontrols.reg <- ivreg(Y ~ 1 + D + X + G | 1 + Z + X + G)
# Output to a LaTeX Tables
stargazer(
    secondstage.reg,
    secondstage_controls.reg,
    secondstage_genecontrols.reg,
    dep.var.caption = "Dependent Variables: Log Annual Income",
    add.lines = list(
        c("Gene Score Controls?", "No", "No", "Yes")),
    #omit = "gene",
    omit.stat = c("LL", "ser", "aic", "adj.rsq"),
    star.cutoffs = NA,
    header = FALSE, float = FALSE, no.space = TRUE,
    omit.table.layout = "n", notes.append = FALSE,
    type = "text",
    out = file.path(tables.folder, "secondstage-reg.tex"))

# Put together all relevant results, one on clean table.

# Output to a LaTeX Tables
stargazer(
    firststage_genecontrols.reg,
    reducedform_genecontrols.reg,
    secondstage_genecontrols.reg,
    #dep.var.caption = "Dependent Variables: Log Annual Income",
    add.lines = list(
        c("Demographic Controls?", "Yes", "Yes", "Yes"),
        c("Gene Score Controls?", "Yes", "Yes", "Yes")),
    omit = "gene|X",
    omit.stat = c("LL", "ser", "aic", "adj.rsq"),
    star.cutoffs = NA,
    header = FALSE, float = FALSE, no.space = TRUE,
    omit.table.layout = "n", notes.append = FALSE,
    type = "text",
    out = file.path(tables.folder, "compiled-reg.tex"))


################################################################################
## New CCI weighting approach.

# Import the pre-defined CCI estimator function, cci.est()
source(file.path("..", "simulations", "cci-define.R"))

# Apply the estimator on the HRS data.
cci.weights <- cci_weights.est(Y, D, Z, cbind(X, G), count.trees = 2000)
cci.est <- cci_point.est(Y, D, Z, cbind(X, G), cci.weights)

# SHow CCI point estimates.
print(c("Total Effect:", cci.est["total_effect"]))
print(c("LAME (for compliers):", cci.est["LAME"]))
print(c("AME:", cci.est["mechanism_effect"]))
print(c("ADE:", cci.est["direct_effect"]))
print(c("Av Indirect E:", cci.est["indirect_effect"]))
print(c("Naive IV Est:", cci.est["iv_est"]))
print(c("IV Bias:", cci.est["iv_bias"]))
# Show the partially identified bounds
print(c("AME bounds:",
    cci.est["mechanism_lower"], cci.est["mechanism_effect"]))

# calculate the SEs by the bootstrap.
cci_bootstrap.data <- cci_bootstrap.est(Y, D, Z, X,
    count.trees = 2000, boot.count = 1000)
tibble(cci_bootstrap.data)

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
estimates.plot <- cci_bootstrap.data %>%
    dplyr::select(#ols_est,
        iv_est, LAME, mechanism_effect) %>%
    pivot_longer(everything(),
        names_to = "parameter", values_to = "estimate") %>%
    ggplot(aes(x = estimate, y = parameter,
        fill = parameter, colour = parameter)) +
    geom_density_ridges2(alpha = 0.5) +
    geom_vline(data = mean_estimate.data,
        aes(xintercept = mean_estimate, colour = parameter),
        linetype = "dashed") +
    theme_bw() +
    scale_x_continuous(name = "Estimate",
        limits = c(0, 1.25),
        breaks = seq(0, 1.25, by = 0.1)) +
    scale_y_discrete(name = "",
        limits = c("mechanism_effect", "LAME", "iv_est"), # "ols_est"),
        breaks = c("iv_est", "LAME", "mechanism_effect"), # "ols_est"),
        labels = c("Naive IV", "CCI (LAME)", "CCI (AME)")) + # "OLS"),
    theme(
        #plot.title = element_text(size = rel(1)),
        plot.margin = unit(c(0.5, 0, 0, 0), "mm"),
        legend.position = "none",
        legend.margin = margin(t = -10))
# Save this plot
ggsave(file.path(figures.folder, "estimates-plot.png"),
    plot = estimates.plot,
    units = "cm", dpi = 300, width = fig.width, height = fig.height)

# Save the SEs
sink(file = file.path(tables.folder, "hrs-cci-ses.txt"))
# Show the boot SEs for the estimates of relevant parameters
print("OLS SEs and 95% CI:")
cci_bootstrap.data %>% pull(ols_est) %>% sd() %>% print()
cci_bootstrap.data %>% pull(ols_est) %>% quantile(probs = c(0.025, 0.975)) %>% print()
print("IV SEs and 95% CI:")
cci_bootstrap.data %>% pull(iv_est) %>% sd() %>% print()
cci_bootstrap.data %>% pull(iv_est) %>% quantile(probs = c(0.025, 0.975)) %>% print()
print("CCI (LAME) SEs and 95% CI:")
cci_bootstrap.data %>% pull(LAME) %>% sd() %>% print()
cci_bootstrap.data %>% pull(LAME) %>% quantile(probs = c(0.025, 0.975)) %>% print()
print("CCI (AME) SEs and 95% CI:")
cci_bootstrap.data %>% pull(mechanism_effect) %>% sd() %>% print()
cci_bootstrap.data %>% pull(mechanism_effect) %>% quantile(probs = c(0.025, 0.975)) %>% print()
# Finish file writing.
sink(file = NULL)


################################################################################
## CCI Conditional estimates.

cci_parentcollege.est <- cci_point.est(Y, D, Z, cbind(X, G), cci.weights,
    subset = (train_hrs_panel.data$parent_edyears >= 16))
cci_noparentcollege.est <- cci_point.est(Y, D, Z, cbind(X, G), cci.weights,
    subset = (train_hrs_panel.data$parent_edyears < 16))
cci_female.est <- cci_point.est(Y, D, Z, cbind(X, G), cci.weights,
    subset = (train_hrs_panel.data$gender_female == 1))
cci_male.est <- cci_point.est(Y, D, Z, cbind(X, G), cci.weights,
    subset = (train_hrs_panel.data$gender_female == 0))
cci_1990s.est <- cci_point.est(Y, D, Z, cbind(X, G), cci.weights,
    subset = (train_hrs_panel.data$survey_year %in% 1990:1999))
cci_2000s.est <- cci_point.est(Y, D, Z, cbind(X, G), cci.weights,
    subset = (train_hrs_panel.data$survey_year %in% 2000:2009))
cci_2010s.est <- cci_point.est(Y, D, Z, cbind(X, G), cci.weights,
    subset = (train_hrs_panel.data$survey_year %in% 2010:2020))

# Store in a dataframe
subset.list <- c("All",
    "Parents College",
    "Parents No College",
    "Female",
    "Male",
    "1990s",
    "2000s",
    "2010s")
est.list <- c(
    cci.est["mechanism_effect"],
        cci.est["direct_effect"],
        cci.est["indirect_effect"],
    cci_parentcollege.est["mechanism_effect"],
        cci_parentcollege.est["direct_effect"],
        cci_parentcollege.est["indirect_effect"],
    cci_noparentcollege.est["mechanism_effect"],
        cci_noparentcollege.est["direct_effect"],
        cci_noparentcollege.est["indirect_effect"],
    cci_female.est["mechanism_effect"],
        cci_female.est["direct_effect"],
        cci_female.est["indirect_effect"],
    cci_male.est["mechanism_effect"],
        cci_male.est["direct_effect"],
        cci_male.est["indirect_effect"],
    cci_1990s.est["mechanism_effect"],
        cci_1990s.est["direct_effect"],
        cci_1990s.est["indirect_effect"],
    cci_2000s.est["mechanism_effect"],
        cci_2000s.est["direct_effect"],
        cci_2000s.est["indirect_effect"],
    cci_2010s.est["mechanism_effect"],
        cci_2010s.est["direct_effect"],
        cci_2010s.est["indirect_effect"])
type.list <- c("Mechanism", "Direct", "Indirect") %>% rep(length(subset.list))

# Put into a dataframe to plot.
conditionalest.data <- data.frame(
    subset = rep(subset.list, each = 3),
    estimate = est.list,
    est_type = type.list) %>%
    tibble()
# Plot a bar chart of the estimates
conditionalest.plot <- conditionalest.data %>%
    ggplot(aes(x = subset, y = estimate, fill = est_type)) +
    geom_col(stat = "identity", position = position_dodge()) +
    coord_flip() +
    theme_bw() +
    scale_y_continuous(name =
        "Estimate",
        #limits = c(-0.025, 0.3),
        breaks = seq(-1, 1, by = 0.1)) +
    scale_x_discrete(name = "", limits = rev(subset.list)) +
    #ggtitle("Conditional Estimates") +
    theme(plot.title = element_text(size = rel(1)),
        plot.margin = unit(c(0.5, 0, 0, 0), "mm"),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.margin = margin(t = -10))
# Save this plot
ggsave(file.path(figures.folder, "conditionalest-plot.png"),
    plot = conditionalest.plot,
    units = "cm", width = fig.width, height = fig.height)


################################################################################
## Compare estimates using causal medation package, Imai Keele Yamamoto (2010)

# Simplify the dataframe
med_hrs_panel.data <- train_hrs_panel.data %>%
    mutate(
        genescore_educ_euro = as.integer(genescore_educ_euro >= 0),
        degree_reported = as.integer(degree_reported >= 5),
        indiv_income_real = log(indiv_income_real))
# Define the first stage / medation model
first_stage.reg <- med_hrs_panel.data %>%
    dplyr::select(- indiv_income_real) %>%
    lm(degree_reported ~ 1 + genescore_educ_euro + ., data = .)
# Define the second stage
second_stage.reg <- lm(indiv_income_real ~
    1 + genescore_educ_euro + degree_reported + ., data = med_hrs_panel.data)
# Estimate the mechanism model.
mechanism.reg <- mediate(first_stage.reg, second_stage.reg,
    treat = "genescore_educ_euro", mediator = "degree_reported",
    robustSE = FALSE, sims = 100)
# Show the mechanism estimates -> they are wrong when D is not independent.
print(cci.est["direct_effect"])
summary(mechanism.reg)
print(c(mechanism.reg$d1 / coef(first_stage.reg)["genescore_educ_euro"],
    cci.est["mechanism_effect"]))
