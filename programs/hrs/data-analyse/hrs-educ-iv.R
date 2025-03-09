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
    filter(!is.na(parent_edyears))
    
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
        parent_edyears, #mother_edyears, father_edyears,
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
#train_hrs_panel.data <- train_hrs_panel.data %>%
#    filtyer(12 <= indiv_edyears & indiv_edyears <= 16)

#! TEST: any missings?
c(dim(drop_na(train_hrs_panel.data)), dim(train_hrs_panel.data))

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
        - starts_with("genescore")) %>%
    as.matrix()
G <- train_hrs_panel.data %>%
    dplyr::select(starts_with("genescore"), - genescore_educ_euro) %>%
    as.matrix()
Y <- train_hrs_panel.data %>%
    transmute(Y = log(indiv_income_real)) %>%
    pull(Y)


################################################################################
## Analyse reduced form + first stage: edu gene score -> income + education.

# naive OLS between educ + income.
wagediff.reg <- lm(Y ~ 1 + D)
wagediff_controls.reg <- lm(Y ~ 1 + D + X)
wagediff_genecontrols.reg <- lm(Y ~ 1 + D + X + G)
summary(wagediff.reg)
summary(wagediff_controls.reg)
#summary(wagediff_genecontrols.reg)

# Linear OLS between educ + gene score
reducedform.reg <- lm(Y ~ 1 + Z)
reducedform_controls.reg <- lm(Y ~ 1 + Z + X)
reducedform_genecontrols.reg <- lm(Y ~ 1 + Z + X + G)
summary(reducedform.reg)
summary(reducedform_controls.reg)
#summary(reducedform_genecontrols.reg)

# Output to a LaTeX Tables
stargazer(
    reducedform.reg,
    reducedform_controls.reg,
    reducedform_genecontrols.reg,
    dep.var.caption = "Dependent Variables: Log Earnings",
    add.lines = list(
        c("Gene Score Controls?", "No", "No", "Yes")),
    omit = "gene",
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
#summary(firststage_genecontrols.reg)

# Naive IV
secondstage.reg <- ivreg(Y ~ 1 + D | 1 + Z)
secondstage_controls.reg <- ivreg(Y ~ 1 + D + X | 1 + Z + X)
secondstage_genecontrols.reg <- ivreg(Y ~ 1 + D + X + G | 1 + Z + X + G)
summary(secondstage.reg)
summary(secondstage_controls.reg)
#summary(secondstage_genecontrols.reg)

# Sequential Selection-on-observables (equivalent to causal mediation)
soo.reg <- lm(Y ~ 1 + D + Z + D:Z)
soo_controls.reg <- lm(Y ~ 1 + D + Z +  + D:Z + X)
soo_genecontrols.reg <- lm(Y ~ 1 + D + Z + D:Z + X + G)
summary(soo.reg)
summary(soo_controls.reg)
#summary(soo_genecontrols.reg)

# Output to a LaTeX Tables
stargazer(
    firststage.reg,
    firststage_controls.reg,
    firststage_genecontrols.reg,
    dep.var.caption = "Dependent Variables: Attended university",
    add.lines = list(
        c("Gene Score Controls?", "No", "No", "Yes")),
    omit = "gene",
    omit.stat = c("LL", "ser", "aic", "adj.rsq"),
    star.cutoffs = NA,
    header = FALSE, float = FALSE, no.space = TRUE,
    omit.table.layout = "n", notes.append = FALSE,
    type = "text",
    out = file.path(tables.folder, "firststage-reg.tex"))

# Output to a LaTeX Tables, for the IV model.
stargazer(
    secondstage.reg,
    secondstage_controls.reg,
    secondstage_genecontrols.reg,
    dep.var.caption = "Dependent Variables: Log Annual Income",
    add.lines = list(
        c("Gene Score Controls?", "No", "No", "Yes")),
    omit = "gene",
    omit.stat = c("LL", "ser", "aic", "adj.rsq"),
    star.cutoffs = NA,
    header = FALSE, float = FALSE, no.space = TRUE,
    omit.table.layout = "n", notes.append = FALSE,
    type = "text",
    out = file.path(tables.folder, "secondstage-reg.tex"))

## Put together all relevant results, one clean table.
# Output to a LaTeX Tables
stargazer(
    firststage_genecontrols.reg,
    reducedform_genecontrols.reg,
    secondstage_genecontrols.reg,
    soo_genecontrols.reg,
    add.lines = list(
        c("Demographic Controls?", "Yes", "Yes", "Yes", "Yes"),
        c("Gene Score Controls?", "Yes", "Yes", "Yes", "Yes")),
    omit = "gene|X",
    omit.stat = c("LL", "ser", "aic", "adj.rsq", "f"),
    star.cutoffs = NA,
    header = FALSE, float = FALSE, no.space = TRUE,
    omit.table.layout = "n", notes.append = FALSE,
    type = "text",
    out = file.path(tables.folder, "compiled-reg.tex"))


################################################################################
## Estimate by causal mediation, with a seuqence of controls.

# Create an empty data frame to store estimates.
totaleffect.est <- c()
totaleffect.lower <- c()
totaleffect.upper <- c()
directeffect.est <- c()
directeffect.lower <- c()
directeffect.upper <- c()
indirecteffect.est <- c()
indirecteffect.lower <- c()
indirecteffect.upper <- c()
percenteffect.est <- c()
percenteffect.lower <- c()
percenteffect.upper <- c()
XG <- cbind(X, G)
X.names <- c(
    "No controls",
    "+ Year",
    "+ Survey Count",
    "+ House count",
    "+ Female",
    "+ Age",
    "+ Parent Edyears",
    "+ Cgnition Score",
    "+ BMI Score",
    "+ Waist Score",
    "+ Height Score",
    "+ Schizophernia Score",
    "+ Alcoholism Score",
    "+ Alzeihmers Score",
    "+ Neurotic Score",
    "+ Wellbeing Score",
    "+ Heart CVD Score",
    "+ Heart MI Score",
    "+ Bipolar Score",
    "+ ADHD Score",
    "+ Cross-disorder Score",
    "+ Extraversion Score",
    "+ Longevity Score",
    "+ Anti-Social Score",
    "+ OCD Score",
    "+ Depression Score",
    "+ PTSD Score",
    "+ Anxiety Factor Score",
    "+ Anxiety Control Score",
    "+ Agesmoking Score")
empty.names <- as.character(0:dim(XG)[2])
# First stage: No controls
# Calculate the causal mediation estimates.
firststage.reg <- lm(D ~ 1 + Z)
soo.reg <- lm(Y ~ 1 + D + Z)
mediation.reg <- mediate(firststage.reg, soo.reg,
    treat = "Z", mediator = "D",
    robustSE = TRUE, sims = 500)
# Extract the total effect estimates.
totaleffect.est <- c(totaleffect.est, summary(mediation.reg)$tau.coef)
totaleffect.lower <- c(totaleffect.lower, as.numeric(summary(mediation.reg)$tau.ci[1]))
totaleffect.upper <- c(totaleffect.upper, as.numeric(summary(mediation.reg)$tau.ci[2]))
# Extract the direct effect estimates.
directeffect.est <- c(directeffect.est, summary(mediation.reg)$z0)
directeffect.lower <- c(directeffect.lower, as.numeric(summary(mediation.reg)$z0.ci[1]))
directeffect.upper <- c(directeffect.upper, as.numeric(summary(mediation.reg)$z0.ci[2]))
# Extract the indirect effect estimates.
indirecteffect.est <- c(indirecteffect.est, summary(mediation.reg)$d0)
indirecteffect.lower <- c(indirecteffect.lower, as.numeric(summary(mediation.reg)$d0.ci[1]))
indirecteffect.upper <- c(indirecteffect.upper, as.numeric(summary(mediation.reg)$d0.ci[2]))
# Extract the percent of total effect mediated estimates.
percenteffect.est <- c(percenteffect.est, summary(mediation.reg)$n.avg)
percenteffect.lower <- c(percenteffect.lower, as.numeric(summary(mediation.reg)$n.avg.ci[1]))
percenteffect.upper <- c(percenteffect.upper, as.numeric(summary(mediation.reg)$n.avg.ci[2]))

# Loop through, adding more controls.
for (i in seq_len(dim(XG)[2])){
    print(i)
    # Calculate the causal mediation estimates.
    firststage.reg <- lm(D ~ 1 + Z + XG[, 1:i])
    soo.reg <- lm(Y ~ 1 + D + Z + XG[, 1:i])
    mediation.reg <- mediate(firststage.reg, soo.reg,
        treat = "Z", mediator = "D",
        robustSE = TRUE, sims = 500)
    # Extract the total effect estimates.
    totaleffect.est <- c(totaleffect.est, summary(mediation.reg)$tau.coef)
    totaleffect.lower <- c(totaleffect.lower, as.numeric(summary(mediation.reg)$tau.ci[1]))
    totaleffect.upper <- c(totaleffect.upper, as.numeric(summary(mediation.reg)$tau.ci[2]))
    # Extract the direct effect estimates.
    directeffect.est <- c(directeffect.est, summary(mediation.reg)$z0)
    directeffect.lower <- c(directeffect.lower, as.numeric(summary(mediation.reg)$z0.ci[1]))
    directeffect.upper <- c(directeffect.upper, as.numeric(summary(mediation.reg)$z0.ci[2]))
    # Extract the indirect effect estimates.
    indirecteffect.est <- c(indirecteffect.est, summary(mediation.reg)$d0)
    indirecteffect.lower <- c(indirecteffect.lower, as.numeric(summary(mediation.reg)$d0.ci[1]))
    indirecteffect.upper <- c(indirecteffect.upper, as.numeric(summary(mediation.reg)$d0.ci[2]))
    # Extract the percent of total effect mediated estimates.
    percenteffect.est <- c(percenteffect.est, summary(mediation.reg)$n.avg)
    percenteffect.lower <- c(percenteffect.lower, as.numeric(summary(mediation.reg)$n.avg.ci[1]))
    percenteffect.upper <- c(percenteffect.upper, as.numeric(summary(mediation.reg)$n.avg.ci[2]))
}

# Store in a datframe
mediation.data <- data.frame(
    X.names = X.names,
    empty.names = empty.names,
    totaleffect.est = totaleffect.est,
    totaleffect.lower = totaleffect.lower,
    totaleffect.upper = totaleffect.upper,
    directeffect.est = directeffect.est,
    directeffect.lower = directeffect.lower,
    directeffect.upper = directeffect.upper,
    indirecteffect.est = indirecteffect.est,
    indirecteffect.lower = indirecteffect.lower,
    indirecteffect.upper = indirecteffect.upper,
    percenteffect.est = percenteffect.est,
    percenteffect.lower = percenteffect.lower,
    percenteffect.upper = percenteffect.upper)

# Plot the total effects.
totaleffect.plot <- mediation.data %>%
    ggplot(aes(x = totaleffect.est, y = X.names)) +
    geom_point(colour = colour.list[1]) +
    geom_errorbar(aes(xmin = totaleffect.lower, xmax = totaleffect.upper,
        width = 0.2), colour = colour.list[1]) +
    theme_bw() +
    theme(legend.position = "none") +
    scale_x_continuous(expand = c(0.01, 0),
        name = "Total Effect Estimate",
        limits = c(0, 0.12),
        breaks = seq(0, 0.15, by = 0.02)) +
    scale_y_discrete(name = "", limits = rev(X.names)) +
    theme(plot.title = element_text(size = rel(1)),
        plot.title.position = "plot",
        plot.margin = unit(c(0.5, 3, 0, 0), "mm"))
# Save this plot
ggsave(file.path(figures.folder, "totaleffect-sequential.png"),
    plot = totaleffect.plot,
    units = "cm", width = 1.25 * fig.width, height = fig.height)
# Plot the indirect effects.
indirecteffect.plot <- mediation.data %>%
    ggplot(aes(x = indirecteffect.est, y = X.names)) +
    geom_point(colour = colour.list[2]) +
    geom_errorbar(aes(xmin = indirecteffect.lower, xmax = indirecteffect.upper,
        width = 0.2), colour = colour.list[2]) +
    theme_bw() +
    theme(legend.position = "none") +
    scale_x_continuous(expand = c(0.01, 0),
        name = "Indirect Effect Estimate",
        limits = c(0, 0.12),
        breaks = seq(0, 0.15, by = 0.02)) +
    scale_y_discrete(name = "", limits = rev(X.names)) +
    theme(plot.title = element_text(size = rel(1)),
        plot.title.position = "plot",
        plot.margin = unit(c(0.5, 3, 0, 0), "mm"))
# Save this plot
ggsave(file.path(figures.folder, "indirecteffect-sequential.png"),
    plot = indirecteffect.plot,
    units = "cm", width = 0.75 * fig.width, height = 1.25 * fig.height)
# Plot the direct effects.
directeffect.plot <- mediation.data %>%
    ggplot(aes(x = directeffect.est, y = empty.names)) +
    geom_point(colour = colour.list[3]) +
    geom_errorbar(aes(xmin = directeffect.lower, xmax = directeffect.upper,
        width = 0.2), colour = colour.list[3]) +
    theme_bw() +
    theme(legend.position = "none") +
    scale_x_continuous(expand = c(0.01, 0),
        name = "Direct Effect Estimate",
        limits = c(0, 0.12),
        breaks = seq(0, 0.15, by = 0.02)) +
    scale_y_discrete(name = "", limits = rev(empty.names)) +
    theme(plot.title = element_text(size = rel(1)),
        plot.title.position = "plot",
        plot.margin = unit(c(0.5, 3, 0, 0), "mm"))
# Save this plot
ggsave(file.path(figures.folder, "directeffect-sequential.png"),
    plot = directeffect.plot,
    units = "cm", width = 0.5 * fig.width, height = 1.25 * fig.height)
# Plot the direct effects.
percenteffect.plot <- mediation.data %>%
    ggplot(aes(x = percenteffect.est, y = X.names)) +
    geom_point(colour = colour.list[1]) +
    geom_errorbar(aes(xmin = percenteffect.lower, xmax = percenteffect.upper,
        width = 0.2), colour = colour.list[1]) +
    theme_bw() +
    theme(legend.position = "none") +
    scale_x_continuous(expand = c(0.01, 0),
        name = "Percent Effect Mediated",
        limits = c(0, 1),
        breaks = seq(0, 1, by = 0.1)) +
    scale_y_discrete(name = "", limits = rev(X.names)) +
    theme(plot.title = element_text(size = rel(1)),
        plot.title.position = "plot",
        plot.margin = unit(c(0.5, 3, 0, 0), "mm"))
# Save this plot
ggsave(file.path(figures.folder, "percenteffect-sequential.png"),
    plot = percenteffect.plot,
    units = "cm", width = 0.75 * fig.width, height = 1.25 * fig.height)
