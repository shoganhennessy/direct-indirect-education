#!/usr/bin/R
## Senan Hogan-Hennessy, 24 March 2025
## UKB data -> First-stage analysis of PGI causal design (done with Ed PGIs).
set.seed(47)
print(Sys.time())

# Functions for data manipulation and visualisation
library(tidyverse)
# Functions for tables into TeX
library(xtable)
# Library for equations in plots
library(latex2exp)
# Library to sumamrise models in plots
library(modelsummary)
# The standard, linear, IV estimator package.
library(ivreg)
# Define number of digits in tables and graphs
digits.no <- 3
# Size for figures
fig.height <- 10
fig.width <- 1.5 * fig.height
presentation.width <- (3 / 2) * fig.width
presentation.height <- presentation.width
# List of 3 default colours.
colour.list <- c(
    "#1f77b4", # Blue
    "#2ca02c", # Green
    "#d62728") # Red
# Define folder path for the general data.
data.folder <- file.path("..", "..", "..", "data", "ukb-restricted", "cleaned")
# Graphics output folder
presentation.folder <- file.path("..", "..", "..", "text", "presentation-files")
figures.folder <- file.path("..", "..", "..", "text", "sections", "figures")
tables.folder <- file.path("..", "..", "..", "text", "sections", "tables")

# Define a function to get the F stat of a specific variable
fstat.get <- function(input.reg, variable.name){
    # Get the variable specific F stats.
    f.test <- car::linearHypothesis(input.reg, test = "F",
        c(paste0(variable.name, "=0")))
    # Extract the statistic.
    f.stat <- f.test["F"] %>%
        unlist() %>%
        nth(2) %>%
        as.numeric() %>%
        round()
    return(f.stat)
}


################################################################################
## Import relevant UKB data. 

# Load the pre-cleaned UKB panel data.
analysis.data <- data.folder %>%
    file.path("ukb-cleaned-pheno.csv") %>%
    read_csv() %>%
    # Get the sibling imputed analysis sample.
    filter(analysis_sample == 1)
# Show the construction
print(analysis.data)
print(names(analysis.data))


################################################################################
## Code up the mediation results, with correlational estimates of ed returns.


control.formula <- paste0("sex_male + recruitedage + sibling_count + urban +",
    "adhd_pgi + asthma_pgi + bipolar_pgi + bmi_pgi + height_pgi +",
    "schizophrenia_pgi + t2diabetes_pgi")


source("oriv.R")
# Causal first-stage
firststage.reg <- GORIV(edyears ~ 1 +
    sex_male + recruitedage + sibling_count + urban + adhd_pgi + asthma_pgi +
    bipolar_pgi + bmi_pgi + height_pgi + schizophrenia_pgi + t2diabetes_pgi,
    "edpgi_all_imputed_self", "edpgi_exclude_imputed_self",
    IV.list = c(
        "edpgi_all_imputed_random", "edpgi_exclude_imputed_random"),
    control.list = c(
        "edpgi_all_imputed_parental", "edpgi_exclude_imputed_parental"),
    IID = "eid", FID = NULL,
    data = analysis.data)
print(summary(firststage.reg))
names(firststage.reg)
# Correlational second-stage
t_secondstage.list <- GORIV(edyears ~ 1 +
    sex_male + recruitedage + sibling_count + urban + adhd_pgi + asthma_pgi +
    bipolar_pgi + bmi_pgi + height_pgi + schizophrenia_pgi + t2diabetes_pgi,
    "edpgi_all_imputed_self", "edpgi_exclude_imputed_self",
    IV.list = c(
        "edpgi_all_imputed_random", "edpgi_exclude_imputed_random"),
    control.list = c(
        "edpgi_all_imputed_parental", "edpgi_exclude_imputed_parental"),
    IID = "eid", FID = NULL,
    data = analysis.data)
#!Testing with my own function.
source("oriv.R")
analysis.data$log_soc_mean_hourly <- log(analysis.data$soc_mean_hourly)
GORIV_mediate("log_soc_mean_hourly", "edyears", control.formula,
    "edpgi_all_imputed_self", "edpgi_exclude_imputed_self",
    IV.list = c(
        "edpgi_all_imputed_random", "edpgi_exclude_imputed_random"),
    #control.list = c(
    #    "edpgi_all_imputed_parental", "edpgi_exclude_imputed_parental"),
    IID = "eid", #FID = "famid",
    data = analysis.data %>% slice_sample(prop = 1, replace = TRUE))


secondstage.reg <- GORIV(log(soc_mean_hourly) ~ 1 +
    edyears +
    sex_male + recruitedage + sibling_count + urban + adhd_pgi + asthma_pgi +
    bipolar_pgi + bmi_pgi + height_pgi + schizophrenia_pgi + t2diabetes_pgi,
    "edpgi_all_imputed_self", "edpgi_exclude_imputed_self",
    IV.list = c(
        "edpgi_all_imputed_random", "edpgi_exclude_imputed_random"),
    control.list = c(
        "edpgi_all_imputed_parental", "edpgi_exclude_imputed_parental"),
    IID = "eid", FID = NULL,
    data = analysis.data)
print(summary(secondstage.reg))
direct.est <- coeftable(oriv_Z_random.reg)["fit_PGI_MAIN", "Estimate"],






analysis.data %>% pull(edpgi_all_imputed_self) %>% summary() %>% print()
analysis.data %>% pull(edpgi_all_imputed_self) %>% sd() %>% print()
analysis.data %>% pull(edpgi_all_imputed_random) %>% summary() %>% print()
analysis.data %>% pull(edpgi_all_imputed_random) %>% var() %>% print()

# Get the demographic possible confounders.
demographic.data <- analysis.data %>%
    select(
        edpgi_all_imputed_self,
        edpgi_all_imputed_random,
        sex_male,
        recruitedage,
        sibling_count,
        urban,
        adhd_pgi,
        asthma_pgi,
        bipolar_pgi,
        bmi_pgi,
        height_pgi,
        schizophrenia_pgi,
        t2diabetes_pgi)
# Label the confounders.
demographic.labels <- c(
    "sex_male"          = "Male",
    "recruitedage"      = "Age",
    "urban"             = "In city",
    "sibling_count"     = "Sibling count",
    "adhd_pgi"          = "Other PGI: ADHD",
    "asthma_pgi"        = "Other PGI: Asthma",
    "bipolar_pgi"       = "Other PGI: Bipolar",
    "bmi_pgi"           = "Other PGI: BMI",
    "t2diabetes_pgi"    = "Other PGI: Diabetes, Type-2",
    "height_pgi"        = "Other PGI: Height",
    "schizophrenia_pgi" = "Other PGI: Schizophrenia")

# Estimate correlations with the regular Ed PGI.
demographic.reg <- demographic.data %>%
    select(- edpgi_all_imputed_random) %>%
    lm(reformulate(".", response = "edpgi_all_imputed_self"), data = .)
# Estimate the correlation with the random component.
demographic_random.reg <- demographic.data %>%
    select(- edpgi_all_imputed_self) %>%
    lm(reformulate(".", response = "edpgi_all_imputed_random"), data = .)

# Plot the first correlations, with the regular Ed PGI.
demographic.plot <- modelplot(demographic.reg,
        coef_map = rev(demographic.labels),
        coef_omit = "Intercept", colour = "blue", size = 1) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    scale_x_continuous(expand = c(0.005, 0.005),
        name = "Correlation Estimate, and 95% Confidence Interval",
        breaks = seq(-1, 1, by = 0.05)) +
    ggtitle(TeX("Demographic Information")) +
    theme(plot.title = element_text(size = rel(1), hjust = 0),
        plot.title.position = "plot",
        plot.margin = unit(c(0.5, 3, 0, 0), "mm"),
        axis.text.y = element_text(hjust = 0))

# Overlay the random component
combined.plot <- list(
    "Ed PGI" = demographic.reg,"Random component" = demographic_random.reg) %>%
    modelplot(., coef_map = rev(demographic.labels),
        coef_omit = "Intercept", size = 1) +
    scale_color_manual(values = colour.list[c(2, 3)]) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    scale_x_continuous(expand = c(0.005, 0.005),
        name = "Correlation Estimate",
        breaks = seq(-1, 1, by = 0.05)) +
    ggtitle(TeX("Demographic Information")) +
    theme(plot.title = element_text(size = rel(1), hjust = 0),
        plot.title.position = "plot",
        plot.margin = unit(c(0.5, 3, 0, 0), "mm"),
        legend.position = "bottom", #c(0.25, 0.95),
        axis.text.y = element_text(hjust = 0))

# Save this plot
ggsave(file.path(figures.folder, "demographic-correlates.png"),
    plot = combined.plot,
    dpi = 300, units = "cm",
    width = fig.width, height = fig.height)


################################################################################
## Testing out the first-stage, with random component of Ed PGI (both sources).
# For notation, write Z for main Ed PGI using Okbay (2022) weights (23&Me)
# and Ztilde for the supplementary Ed PGI from UKB subsample.

# Random component Z -> Z
Z_Z_random.reg <- lm(edpgi_all_imputed_self ~ 1 +
    sex_male + recruitedage + sibling_count + urban + adhd_pgi + asthma_pgi +
    bipolar_pgi + bmi_pgi + height_pgi + schizophrenia_pgi + t2diabetes_pgi +
    edpgi_all_imputed_random + edpgi_all_imputed_parental,
    data = analysis.data)
Z_Z_random.reg %>% summary() %>% print()

# Random component Ztilde -> Z
Ztilde_Z_random.reg <- lm(edpgi_all_imputed_self ~ 1 +
    sex_male + recruitedage + sibling_count + urban + adhd_pgi + asthma_pgi +
    bipolar_pgi + bmi_pgi + height_pgi + schizophrenia_pgi + t2diabetes_pgi + 
    edpgi_exclude_imputed_random + edpgi_exclude_imputed_parental,
    data = analysis.data)
Ztilde_Z_random.reg %>% summary() %>% print()

# Random component Z -> Ztilde
Z_Ztilde_random.reg <- lm(edpgi_exclude_imputed_self ~ 1 +
    sex_male + recruitedage + sibling_count + urban + adhd_pgi + asthma_pgi +
    bipolar_pgi + bmi_pgi + height_pgi + schizophrenia_pgi + t2diabetes_pgi + 
    edpgi_all_imputed_random + edpgi_all_imputed_parental,
    data = analysis.data)
Z_Ztilde_random.reg %>% summary() %>% print()

# Random component Ztilde -> Ztilde
Ztilde_Ztilde_random.reg <- lm(edpgi_exclude_imputed_self ~ 1 +
    sex_male + recruitedage + sibling_count + urban + adhd_pgi + asthma_pgi +
    bipolar_pgi + bmi_pgi + height_pgi + schizophrenia_pgi + t2diabetes_pgi +
    edpgi_exclude_imputed_random + edpgi_exclude_imputed_parental,
    data = analysis.data)
Ztilde_Ztilde_random.reg %>% summary() %>% print()

## Estimate as an ORIV first-stage (stacked).
source("oriv.R")
# Z main
oriv_Z_random.reg <- GORIV(edpgi_all_imputed_self ~ 1 +
    sex_male + recruitedage + sibling_count + urban + adhd_pgi + asthma_pgi +
    bipolar_pgi + bmi_pgi + height_pgi + schizophrenia_pgi + t2diabetes_pgi,
    "edpgi_all_imputed_random", "edpgi_exclude_imputed_random",
    IV.list = c(
        "edpgi_all_imputed_random", "edpgi_exclude_imputed_random"),
    control.list = c(
        "edpgi_all_imputed_parental", "edpgi_exclude_imputed_parental"),
    IID = "eid",
    data = analysis.data)
oriv_Z_random.reg %>% summary() %>% print()
# Z tilde
oriv_Ztilde_random.reg <- GORIV(edpgi_exclude_imputed_self ~ 1 +
    sex_male + recruitedage + sibling_count + urban + adhd_pgi + asthma_pgi +
    bipolar_pgi + bmi_pgi + height_pgi + schizophrenia_pgi + t2diabetes_pgi,
    "edpgi_all_imputed_random", "edpgi_exclude_imputed_random",
    IV.list = c(
        "edpgi_all_imputed_random", "edpgi_exclude_imputed_random"),
    control.list = c(
        "edpgi_all_imputed_parental", "edpgi_exclude_imputed_parental"),
    IID = "eid",
    data = analysis.data)
oriv_Ztilde_random.reg %>% summary() %>% print()


################################################################################
## LaTeX Table of the genetic first-stage.

# Compose all the point estimates and SEs, then F stat, obs count, R^2.
point_est.data <- data.frame(
    # column 1: Z -> Z
    col1 = c(
        coef(summary(Z_Z_random.reg))["edpgi_all_imputed_random", "Estimate"],
        coef(summary(Z_Z_random.reg))["edpgi_all_imputed_random", "Std. Error"],
        NA, NA, NA, NA,
        NA, # fstat.get(Z_Z_random.reg, "edpgi_all_imputed_random"),
        summary(Z_Z_random.reg)$r.squared,
        NROW(analysis.data)),
    # column 2: Ztilde -> Z
    col2 = c(NA, NA,
        coef(summary(Ztilde_Z_random.reg))["edpgi_exclude_imputed_random", "Estimate"],
        coef(summary(Ztilde_Z_random.reg))["edpgi_exclude_imputed_random", "Std. Error"],
        NA, NA,
        fstat.get(Ztilde_Z_random.reg, "edpgi_exclude_imputed_random"),
        summary(Ztilde_Z_random.reg)$r.squared,
        NROW(analysis.data)),
    # column 3: obviously related Z, Ztilde -> Z
    col3 = c(NA, NA, NA, NA,
        coeftable(oriv_Z_random.reg)["fit_PGI_MAIN", "Estimate"],
        coeftable(oriv_Z_random.reg)["fit_PGI_MAIN", "Std. Error"],
        fstat.get(oriv_Z_random.reg, "fit_PGI_MAIN"),
        NA,
        NROW(analysis.data)),
    # Column 4: Z -> Ztilde
    col4 = c(
        coef(summary(Z_Ztilde_random.reg))["edpgi_all_imputed_random", "Estimate"],
        coef(summary(Z_Ztilde_random.reg))["edpgi_all_imputed_random", "Std. Error"],
        NA, NA, NA, NA,
        fstat.get(Z_Ztilde_random.reg, "edpgi_all_imputed_random"),
        summary(Z_Ztilde_random.reg)$r.squared,
        NROW(analysis.data)),
    # column 5: Ztilde -> Ztilde
    col5 = c(NA, NA,
        coef(summary(Ztilde_Ztilde_random.reg))["edpgi_exclude_imputed_random", "Estimate"],
        coef(summary(Ztilde_Ztilde_random.reg))["edpgi_exclude_imputed_random", "Std. Error"],
        NA, NA,
        NA, # fstat.get(Ztilde_Ztilde_random.reg, "edpgi_exclude_imputed_random"),
        summary(Ztilde_Ztilde_random.reg)$r.squared,
        NROW(analysis.data)),
    # column 6: obviously related Z, Ztilde -> Z
    col6 = c(NA, NA, NA, NA,
        coeftable(oriv_Ztilde_random.reg)["fit_PGI_MAIN", "Estimate"],
        coeftable(oriv_Ztilde_random.reg)["fit_PGI_MAIN", "Std. Error"],
        fstat.get(oriv_Ztilde_random.reg, "fit_PGI_MAIN"),
        NA,
        NROW(analysis.data))
    )

# Show what you get.
print(point_est.data)

# Round and remove NAs.
point_est.data <- round(point_est.data, digits.no)
point_est.data <- data.frame(lapply(point_est.data, as.character), stringsAsFactors = FALSE) %>%
    replace_na(list(
        col1 = " ",
        col2 = " ",
        col3 = " ",
        col4 = " ",
        col5 = " ",
        col6 = " "))

# Bracket the SEs
for (col.index in 1:3){
    point_est.data[2 * col.index, col.index] <- 
        point_est.data[2 * col.index, col.index]  %>%
        as.numeric() %>%
        signif(digits.no) %>%
        format(scientific = FALSE) %>%
        paste0("(", ., ")")
    point_est.data[2 * col.index, col.index + 3] <-
        point_est.data[2 * col.index, col.index + 3]  %>%
        as.numeric() %>%
        signif(digits.no) %>%
        format(scientific = FALSE) %>%
        paste0("(", ., ")")
}
# Add on variable names (first column)
point_est.data$names = c(
    r"(Primary Ed PGI, $Z_i^\text{\textcolor{red}{Random}}$)", " ",
    r"(Secondary Ed PGI, $Z_i^\text{\textcolor{red}{Random}}$)", " ",
    "Combined, ORIV", " ",
    r"(\\[-1.8ex]\hline \\[-1.8ex] $F$-statistics)",
    "$R^2$",
    "Observations")
point_est.data <- point_est.data[, c(7, 1:6)]

# Show the LaTeX table
point_est.data %>%
    xtable() %>%
    print(
        digits = digits.no,
        sanitize.colnames.function = identity,
        sanitize.text.function = identity,
        NA.string = " ",
        include.colnames = FALSE,
        include.rownames = FALSE,
        only.contents = TRUE,
        hline.after = NULL,
        format.args = list(big.mark = ","),
        file = file.path(tables.folder, "genetic-firststage.tex"))


################################################################################
## Genetic effects

## OLS versus random component.
lm(edyears ~ 1 + edpgi_all_imputed_self, data = analysis.data) %>% summary()
lm(edyears ~ 1 + edpgi_all_imputed_self + edpgi_all_imputed_parental, data = analysis.data) %>% summary()
ivreg(edyears ~ 1 + edpgi_all_imputed_self + edpgi_all_imputed_parental |
    edpgi_all_imputed_random  + edpgi_all_imputed_parental, data = analysis.data) %>% summary()

ivreg(edyears ~ 1 + edpgi_all_imputed_self + edpgi_all_imputed_parental |
    edpgi_exclude_imputed_random  + edpgi_exclude_imputed_parental, data = analysis.data) %>% summary()


ivreg(as.integer(edyears >= 18) ~ 1 + edpgi_all_imputed_self + edpgi_all_imputed_parental |
    edpgi_all_imputed_random + edpgi_all_imputed_parental, data = analysis.data) %>% summary()
#! About twice as large as literature estimates, when using diferent weights for imputed random/parental component.


# Family FEs?
fam.data <- analysis.data %>%
    group_by(famid) %>%
    mutate(fam_count = 1) %>%
    mutate(fam_count = sum(fam_count)) %>%
    ungroup() %>%
    filter(fam_count > 1) %>%
    select(-fam_count)

# Base line: which Ed PGI is most predictive (within families) of Ed years?
feols(edyears ~ 1 + edpgi_all_imputed_self, data = analysis.data) %>% summary()
feols(edyears ~ 1 + edpgi_all_imputed_self | famid,
    #famid | edpgi_all_imputed_self ~ edpgi_all_imputed_random,
    data = fam.data) %>% summary()
# Ed PGI (all raw)              gives 1.09, 0.54
# Phased Ed PGI (all imputed)   gives 1.1, 0.54
# Unphased Ed PGI (all imputed) gives ?, ?
# Ed PGI (exclude raw)          gives 0.186, 0.111
# Ed PGI (exclude imputed)      gives 0.19, 0.116
#! Conclusion: use the unphased Ed PGI.
fam.data


lm(edpgi_all_imputed_self ~ 1 + edpgi_exclude_imputed_random, data = analysis.data) %>% summary()



ivreg(as.integer(edyears >= 18) ~ 1 + edpgi_all_imputed_self + edpgi_all_imputed_parental | edpgi_all_imputed_random + edpgi_all_imputed_parental, data = analysis.data) %>% summary()

ivreg(as.integer(edyears >= 18) ~ 1 + edpgi_all_imputed_self | edpgi_all_imputed_random , data = analysis.data) %>% summary()


# The gain from ORIV (correlation)
summary(lm(edyears ~ 1 + edpgi_all_imputed_self, data = analysis.data))
summary(GORIV(edyears ~ 1,
    "edpgi_all_imputed_self", "edpgi_exclude_imputed_self",
    data = analysis.data, IID = "eid"))

# The gain from ORIV (causal)
source("oriv.R")
summary(lm(edyears ~ 1 + edpgi_all_imputed_self + edpgi_all_imputed_parental,
    data = analysis.data))
summary(GORIV(edyears ~ 1,
    "edpgi_all_imputed_self", "edpgi_exclude_imputed_self",
    IV.list = c(
        "edpgi_all_imputed_random", "edpgi_exclude_imputed_random"),
    control.list =
        c("edpgi_all_imputed_parental", "edpgi_exclude_imputed_parental"),
    data = analysis.data, IID = "eid"))

# For higher ed attendance
summary(lm(I(edyears >= 18) ~ 1 + edpgi_all_imputed_self + edpgi_all_imputed_parental,
    data = analysis.data))
summary(GORIV(I(edyears >= 18) ~ 1,
    "edpgi_all_imputed_self", "edpgi_exclude_imputed_self",
    IV.list = c(
        "edpgi_all_imputed_random", "edpgi_exclude_imputed_random"),
    control.list =
        c("edpgi_all_imputed_parental", "edpgi_exclude_imputed_parental"),
    data = analysis.data, IID = "eid"))




# Causal, with + with/out ORIV measurement error adjustment.
summary(lm(edyears ~ 1 + edpgi_all_imputed_random, data = analysis.data))
summary(ivreg(edyears ~ 1 + edpgi_all_imputed_self | edpgi_all_imputed_random , data = analysis.data))
summary(GORIV(edyears ~ 1,
    "edpgi_all_imputed_random", "edpgi_exclude_imputed_random",
    IID = "eid", FID = "famid",
    data = analysis.data))
# Causal, with + with/out ORIV measurement error adjustment.
summary(lm(log(soc_mean_hourly) ~ 1 + edpgi_all_imputed_random, data = analysis.data))
summary(GORIV(log(soc_mean_hourly) ~ 1,
    "edpgi_all_imputed_random", "edpgi_exclude_imputed_random",
    data = analysis.data, IID = "eid", FID = "famid"))


################################################################################
## Quasi second-stage, Ed PGI -> Ed Years

# OLS (not causal)
analysis.data %>%
    lm(edyears ~ 1 + edpgi_self, data = .) %>%
    summary() %>%
    print()

# Show the correlation between the random value, and Ed years
analysis.data %>%
    lm(edyears ~ 1 + edpgi_random, data = .) %>%
    summary() %>%
    print()

# Show the correlation between the self + parents value, and Ed years
analysis.data %>%
    lm(edyears ~ 1 + edpgi_self + edpgi_parents, data = .) %>%
    summary() %>%
    print()

# Try it with the random part as an instrument.
analysis.data %>%
    ivreg::ivreg(edyears ~ 1 + edpgi_self | 1 + edpgi_random,
        data = .) %>%
    summary() %>%
    print()

# Show the mean Ed PGI random component across the parent dist
analysis.data$parental_quantile <-
    ecdf(analysis.data$edpgi_parents)(analysis.data$edpgi_parents)
analysis.data$parental_quantile <- factor(
    round(10 * analysis.data$parental_quantile))
# Compare this to Houmark+ (2024) Figure 1.
analysis.data %>%
    group_by(parental_quantile) %>%
    summarise(
        edpgi_random_mean   = mean(edpgi_random),
        edpgi_random_sd     = sd(edpgi_random),
        edpgi_self_mean     = mean(edpgi_self),
        edpgi_self_sd       = sd(edpgi_self),
        edpgi_parental_mean = mean(edpgi_parents),
        edpgi_parental_sd   = sd(edpgi_parents)) %>%
    View()
mean(analysis.data$edpgi_self > analysis.data$edpgi_parents)
# Ensure this does not vary across the parental distribution.
analysis.data %>%
    lm(edpgi_self ~ 1 + edpgi_random * parental_quantile, data = .) %>%
    summary() %>%
    print()
#  TEST: does this association hold true among people for whom we observe one or both parents? -> Yes.
analysis.data %>%
    filter(father_present + mother_present > 0) %>%
    lm(edpgi_random ~ 1 + 
        sex_male * (edpgi_father + edpgi_mother), data = .) %>%
    summary() %>%
    print()
# Compare to the Ed PGI in raw form.
analysis.data %>%
    filter(father_present + mother_present > 0) %>%
    lm(edpgi_self ~ 1 +
        sex_male * (edpgi_father + edpgi_mother), data = .) %>%
    summary() %>%
    print()
