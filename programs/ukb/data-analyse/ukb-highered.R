#!/usr/bin/R
## Senan Hogan-Hennessy, 30 July 2025
## UKB data, Effect of Rosla -> Edyears. 
set.seed(47)
print(Sys.time())

# Functions for data manipulation and visualisation
library(tidyverse)
# Functions for tables into TeX
library(xtable)
# Library for equations in plots
library(latex2exp)
# Library to sumamrise models in plots
library(gtsummary)
library(modelsummary)
# The standard, linear, IV estimator package.
library(fixest)
library(ivreg)
# Define number of digits in tables and graphs
digits.no <- 3
# Size for figures
fig.height <- 8.5
fig.width <- 1.25 * fig.height
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

# Load the pre-cleaned UKB cross data.
analysis.data <- data.folder %>%
    file.path("ukb-cleaned-pheno.csv") %>%
    read_csv() %>%
    # Get the sibling imputed analysis sample.
    filter(analysis_sample == 1)


################################################################################
## Higher education returns.

## Log Income: Estimate every model necessary
# 1. Ed years, OLS
ols_edyears_wage.reg <- analysis.data %>%
    feols(log(soc_median_hourly) ~ 1 + edyears | birthyear,
        data = .)
print(summary(ols_edyears_wage.reg))
# 1.5. Ed years, Sibling FE
fe_edyears_wage.reg <- analysis.data %>%
    feols(log(soc_median_hourly) ~ 1 + edyears | birthyear + famid,
        se = "twoway", cluster = c("birthyear", "famid"),
        data = .)
print(summary(fe_edyears_wage.reg))
# 2. Uni degree, OLS
ols_highered_wage.reg <- analysis.data %>%
    filter(12 <= edyears) %>%
    feols(log(soc_median_hourly) ~ 1 + edqual_highered | birthyear,
        data = .)
print(summary(ols_highered_wage.reg))
# 2.5. Ed years, Sibling FE
fe_highered_wage.reg <- analysis.data %>%
    filter(12 <= edyears) %>%
    feols(log(soc_median_hourly) ~ 1 + edqual_highered | birthyear + famid,
        se = "twoway", cluster = c("birthyear", "famid"),
        data = .)
print(summary(fe_highered_wage.reg))

## Log Income: Estimate every model necessary
# 1. Ed years, OLS
ols_edyears_income.reg <- analysis.data %>%
    feols(log(soc_median_annual) ~ 1 + edyears | birthyear,
        data = .)
print(summary(ols_edyears_income.reg))
# 1.5. Ed years, Sibling FE
fe_edyears_income.reg <- analysis.data %>%
    feols(log(soc_median_annual) ~ 1 + edyears | birthyear + famid,
        se = "twoway", cluster = c("birthyear", "famid"),
        data = .)
print(summary(fe_edyears_income.reg))
# 2. Uni degree, OLS
ols_highered_income.reg <- analysis.data %>%
    filter(12 <= edyears) %>%
    feols(log(soc_median_annual) ~ 1 + edqual_highered | birthyear,
        data = .)
print(summary(ols_highered_income.reg))
# 2.5. Ed years, Sibling FE
fe_highered_income.reg <- analysis.data %>%
    filter(12 <= edyears) %>%
    feols(log(soc_median_annual) ~ 1 + edqual_highered | birthyear + famid,
        se = "twoway", cluster = c("birthyear", "famid"),
        data = .)
print(summary(fe_highered_income.reg))

# Collect into a list.
edyears.list <- list(
    summary(ols_edyears_wage.reg),
    summary(ols_edyears_income.reg),
    summary(fe_edyears_wage.reg),
    summary(fe_edyears_income.reg))
highered.list <- list(
    summary(ols_highered_wage.reg),
    summary(ols_highered_income.reg),
    summary(fe_highered_wage.reg),
    summary(fe_highered_income.reg))
order.list <- c(edyears.list, highered.list)[c(1:2, 5:6, 3:4, 7:8)]

# Keep everything same format.
edyears.point <- sapply(edyears.list, function(x) coeftable(x)["edyears", "Estimate"])
edyears.point <- style_sigfig(edyears.point, digits = digits.no, decimals = digits.no)
edyears.point <- as.vector(rbind(edyears.point, " "))
highered.point <- sapply(highered.list, function(x) coeftable(x)["edqual_highered", "Estimate"])
highered.point <- style_sigfig(highered.point, digits = digits.no, decimals = digits.no)
highered.point <- as.vector(rbind(" ", highered.point))
point.obs <- sapply(order.list,
    function(x) format(x$nobs, big.mark = ",", scientific = FALSE))
# Get SEs
edyears.se <- sapply(edyears.list, function(x) coeftable(x)["edyears", "Std. Error"])
edyears.se <- style_sigfig(edyears.se, digits = digits.no, decimals = digits.no) %>% paste0("(", ., ")")
edyears.se <- as.vector(rbind(edyears.se, " "))
highered.se <- sapply(highered.list, function(x) coeftable(x)["edqual_highered", "Std. Error"])
highered.se <- style_sigfig(highered.se, digits = digits.no, decimals = digits.no) %>% paste0("(", ., ")")
highered.se <- as.vector(rbind(" ", highered.se))
# Get R^2
point.r2 <- sapply(order.list,
    function(x) style_sigfig(x$sq.cor,
        digits = digits.no, decimals = digits.no))
# Label the Sibling FEs
point.fes <- rep(c(" ", "Yes"), 4)

# Build a LaTeX table.
returns.table <- rbind(edyears.point, edyears.se, highered.point, highered.se)
# Add R^2, Sibling FEs, obs count.
returns.table <- returns.table %>% rbind(point.r2)
#returns.table <- returns.table %>% rbind(point.fes)
returns.table <- returns.table %>% rbind(point.obs)

# Give it names.
returns.table <- c("Education years", " ", "University degree", " ", 
        r"(\\[-1.8ex]\hline \\[-1.8ex] $R^2$)",
        "Observations") %>%
    cbind(returns.table)

# Show the LaTeX table
returns.table %>%
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
        file = file.path(tables.folder, "highered-returns.tex"))


################################################################################
## Sibling FEs for mediation effects.



library(mediation)
library(mgcv)

analysis.data$edyears_factor <- factor(
    ifelse(analysis.data$edqual_highered == 1, 3,
        ifelse(analysis.data$edqual_alevels == 1, 2,
            ifelse(analysis.data$edqual_gcses == 1, 1, 0))))

analysis.data$log_soc_median_hourly <- log(analysis.data$soc_median_hourly)

firststage.reg <- lm(edyears ~ 1 + edpgi_all_imputed_self
    + edpgi_all_imputed_parental,
    data = analysis.data)
secondstage.reg <- lm(log_soc_median_hourly ~ 1
    + (edpgi_all_imputed_self + edpgi_all_imputed_parental) * edyears,
    data = analysis.data)

print(summary(secondstage.reg))

mediate.reg <- mediate(
    firststage.reg,  # First-stage
    secondstage.reg, # Second-stage
    treat = "edpgi_all_imputed_self", mediator = "edyears",
    sims = 50)
print(summary(mediate.reg))


analysis.data <- analysis.data %>%
    filter(abs(birthyearmonth - cutoff) <= bandwidth) %>%
    mutate(log_soc_median_hourly = log(soc_median_hourly))
dist <- abs((analysis.data$birthyearmonth - cutoff) / bandwidth)
analysis.data$kernel.wt <- (1 - dist) * (dist <= 1) / bandwidth

#! Bootstrap test.
full_boot.data <- full.data
analysis_boot.data <- analysis.data
#full_boot.data <- sample_frac(full.data, 1, replace = TRUE)
#analysis_boot.data <- sample_frac(analysis.data, 1, replace = TRUE)

source("oriv.R")
# Full sample IV first-stage.
rosla.reg <- full_boot.data %>%
    feols(edyears ~ 1 + rosla +
        poly(rosla_year_below, poly.count) + poly(rosla_year_above, poly.count),
        weights = full.data$kernel.wt,
        data = .)
rosla.reg %>% summary() %>% print()
rosla.reg %>% fstat.get("rosla") %>% print()

# Save the fitted values, for a 2SLS approach.
analysis_boot.data$fit_edyears <- predict(rosla.reg, newdata = analysis_boot.data)

# Mediation total effect.
total.reg <- feols(log_soc_median_hourly ~ 1 +
    edpgi_all_imputed_self + edpgi_all_imputed_parental +
    poly(rosla_year_below, poly.count) + poly(rosla_year_above, poly.count),
    #| famid,
    #weights = analysis_boot.data$kernel.wt,
    data = analysis_boot.data)
print(summary(total.reg))

# Mediation first-stage
firststage.reg <- analysis_boot.data %>%
    #filter(edyears >= 10) %>%
    lm(edqual_alevels ~ 1 +
        edpgi_all_imputed_self + edpgi_all_imputed_parental +
        poly(rosla_year_below, poly.count) + poly(rosla_year_above, poly.count),
        #| famid,
        #weights = analysis_boot.data$kernel.wt,
        data = .)
print(summary(firststage.reg))

# Mediation second-stage (with 2SLS predicted mediation).
secondstage.reg <- analysis_boot.data %>%
    #filter(edyears >= 10) %>%
    lm(log_soc_median_hourly ~ 1 +
        edpgi_all_imputed_self * edqual_alevels +
        edpgi_all_imputed_parental +
        poly(rosla_year_below, poly.count) + poly(rosla_year_above, poly.count),
        #| famid,
        #weights = analysis_boot.data$kernel.wt,
        data = .)
print(summary(secondstage.reg))

library(mediation)
mediate.reg <- mediate(
    firststage.reg,  # First-stage
    secondstage.reg, # Second-stage
    treat = "edpgi_all_imputed_self", mediator = "edqual_alevels",
    robustSE = TRUE, sims = 500)
summary(mediate.reg)


mediator <- "edyears"
firststage.est <- coeftable(firststage.reg)["edpgi_all_imputed_self", "Estimate"]
total.est <- coeftable(total.reg)["edpgi_all_imputed_self", "Estimate"]
direct.est <- coeftable(secondstage.reg)["edpgi_all_imputed_self", "Estimate"]
indirect.est <- coeftable(secondstage.reg)[mediator, "Estimate"]
interaction.est <- coeftable(secondstage.reg)[paste0("edpgi_all_imputed_self:", mediator), "Estimate"]
ade.est <- direct.est + interaction.est * mean(
    analysis_boot.data[["edqual_gcses"]], na.rm = TRUE)
aie.est <- firststage.est * (indirect.est +
    interaction.est * mean(analysis_boot.data[["edpgi_all_imputed_self"]]))    
print(c(total.est, firststage.est, direct.est, indirect.est, interaction.est,
    ade.est, aie.est, ade.est / total.est, aie.est / total.est))



# Outcomes: ADE for edyears is roughly equal to the ATE, even when using the IV.
# The effects are roughly similar when using correlational terms because the IV estimates are similarl.

# The ADE is smaller for correlational estimates of higher ed,
# implying at least 50\% goes through higher education.

# The next thing to do here would be a sequential model of education choice,
# Choose to finish GCSES, A-Levels, and then unviersity.
# Which emits then a structural solution to the indirect effect through GCSES,
# then indirect effect via A-Levels (given GCSES finished),
# then indirect effect via uni degree (given A-Levels finished).



#!Note, this will require doing a split sample RDD
#! Use the entire sample in the first-stage, predicting education from 

#! ALso note that this is not robust to a polynomial different from 2.



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
    data = full.data)
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
    data = full.data)
#!Testing with my own function.
source("oriv.R")
full.data$log_soc_median_hourly <- log(full.data$soc_median_hourly)
GORIV_mediate("log_soc_median_hourly", "edyears", control.formula,
    "edpgi_all_imputed_self", "edpgi_exclude_imputed_self",
    IV.list = c(
        "edpgi_all_imputed_random", "edpgi_exclude_imputed_random"),
    #control.list = c(
    #    "edpgi_all_imputed_parental", "edpgi_exclude_imputed_parental"),
    IID = "eid", #FID = "famid",
    data = full.data %>% slice_sample(prop = 1, replace = TRUE))


secondstage.reg <- GORIV(log(soc_median_hourly) ~ 1 +
    edyears +
    sex_male + recruitedage + sibling_count + urban + adhd_pgi + asthma_pgi +
    bipolar_pgi + bmi_pgi + height_pgi + schizophrenia_pgi + t2diabetes_pgi,
    "edpgi_all_imputed_self", "edpgi_exclude_imputed_self",
    IV.list = c(
        "edpgi_all_imputed_random", "edpgi_exclude_imputed_random"),
    control.list = c(
        "edpgi_all_imputed_parental", "edpgi_exclude_imputed_parental"),
    IID = "eid", FID = NULL,
    data = full.data)
print(summary(secondstage.reg))
direct.est <- coeftable(oriv_Z_random.reg)["fit_PGI_MAIN", "Estimate"],






full.data %>% pull(edpgi_all_imputed_self) %>% summary() %>% print()
full.data %>% pull(edpgi_all_imputed_self) %>% sd() %>% print()
full.data %>% pull(edpgi_all_imputed_random) %>% summary() %>% print()
full.data %>% pull(edpgi_all_imputed_random) %>% var() %>% print()

# Get the demographic possible confounders.
demographic.data <- full.data %>%
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
    data = full.data)
Z_Z_random.reg %>% summary() %>% print()

# Random component Ztilde -> Z
Ztilde_Z_random.reg <- lm(edpgi_all_imputed_self ~ 1 +
    sex_male + recruitedage + sibling_count + urban + adhd_pgi + asthma_pgi +
    bipolar_pgi + bmi_pgi + height_pgi + schizophrenia_pgi + t2diabetes_pgi + 
    edpgi_exclude_imputed_random + edpgi_exclude_imputed_parental,
    data = full.data)
Ztilde_Z_random.reg %>% summary() %>% print()

# Random component Z -> Ztilde
Z_Ztilde_random.reg <- lm(edpgi_exclude_imputed_self ~ 1 +
    sex_male + recruitedage + sibling_count + urban + adhd_pgi + asthma_pgi +
    bipolar_pgi + bmi_pgi + height_pgi + schizophrenia_pgi + t2diabetes_pgi + 
    edpgi_all_imputed_random + edpgi_all_imputed_parental,
    data = full.data)
Z_Ztilde_random.reg %>% summary() %>% print()

# Random component Ztilde -> Ztilde
Ztilde_Ztilde_random.reg <- lm(edpgi_exclude_imputed_self ~ 1 +
    sex_male + recruitedage + sibling_count + urban + adhd_pgi + asthma_pgi +
    bipolar_pgi + bmi_pgi + height_pgi + schizophrenia_pgi + t2diabetes_pgi +
    edpgi_exclude_imputed_random + edpgi_exclude_imputed_parental,
    data = full.data)
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
    data = full.data)
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
    data = full.data)
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
        NROW(full.data)),
    # column 2: Ztilde -> Z
    col2 = c(NA, NA,
        coef(summary(Ztilde_Z_random.reg))["edpgi_exclude_imputed_random", "Estimate"],
        coef(summary(Ztilde_Z_random.reg))["edpgi_exclude_imputed_random", "Std. Error"],
        NA, NA,
        fstat.get(Ztilde_Z_random.reg, "edpgi_exclude_imputed_random"),
        summary(Ztilde_Z_random.reg)$r.squared,
        NROW(full.data)),
    # column 3: obviously related Z, Ztilde -> Z
    col3 = c(NA, NA, NA, NA,
        coeftable(oriv_Z_random.reg)["fit_PGI_MAIN", "Estimate"],
        coeftable(oriv_Z_random.reg)["fit_PGI_MAIN", "Std. Error"],
        fstat.get(oriv_Z_random.reg, "fit_PGI_MAIN"),
        NA,
        NROW(full.data)),
    # Column 4: Z -> Ztilde
    col4 = c(
        coef(summary(Z_Ztilde_random.reg))["edpgi_all_imputed_random", "Estimate"],
        coef(summary(Z_Ztilde_random.reg))["edpgi_all_imputed_random", "Std. Error"],
        NA, NA, NA, NA,
        fstat.get(Z_Ztilde_random.reg, "edpgi_all_imputed_random"),
        summary(Z_Ztilde_random.reg)$r.squared,
        NROW(full.data)),
    # column 5: Ztilde -> Ztilde
    col5 = c(NA, NA,
        coef(summary(Ztilde_Ztilde_random.reg))["edpgi_exclude_imputed_random", "Estimate"],
        coef(summary(Ztilde_Ztilde_random.reg))["edpgi_exclude_imputed_random", "Std. Error"],
        NA, NA,
        NA, # fstat.get(Ztilde_Ztilde_random.reg, "edpgi_exclude_imputed_random"),
        summary(Ztilde_Ztilde_random.reg)$r.squared,
        NROW(full.data)),
    # column 6: obviously related Z, Ztilde -> Z
    col6 = c(NA, NA, NA, NA,
        coeftable(oriv_Ztilde_random.reg)["fit_PGI_MAIN", "Estimate"],
        coeftable(oriv_Ztilde_random.reg)["fit_PGI_MAIN", "Std. Error"],
        fstat.get(oriv_Ztilde_random.reg, "fit_PGI_MAIN"),
        NA,
        NROW(full.data))
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
        style_sigfig(digits = digits.no) %>%
        format(scientific = FALSE) %>%
        paste0("(", ., ")")
    point_est.data[2 * col.index, col.index + 3] <-
        point_est.data[2 * col.index, col.index + 3]  %>%
        as.numeric() %>%
        style_sigfig(digits = digits.no) %>%
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
lm(edyears ~ 1 + edpgi_all_imputed_self, data = full.data) %>% summary()
lm(edyears ~ 1 + edpgi_all_imputed_self + edpgi_all_imputed_parental, data = full.data) %>% summary()
ivreg(edyears ~ 1 + edpgi_all_imputed_self + edpgi_all_imputed_parental |
    edpgi_all_imputed_random  + edpgi_all_imputed_parental, data = full.data) %>% summary()

ivreg(edyears ~ 1 + edpgi_all_imputed_self + edpgi_all_imputed_parental |
    edpgi_exclude_imputed_random  + edpgi_exclude_imputed_parental, data = full.data) %>% summary()


ivreg(as.integer(edyears >= 18) ~ 1 + edpgi_all_imputed_self + edpgi_all_imputed_parental |
    edpgi_all_imputed_random + edpgi_all_imputed_parental, data = full.data) %>% summary()
#! About twice as large as literature estimates, when using diferent weights for imputed random/parental component.


# Family FEs?
fam.data <- full.data %>%
    group_by(famid) %>%
    mutate(fam_count = 1) %>%
    mutate(fam_count = sum(fam_count)) %>%
    ungroup() %>%
    filter(fam_count > 1) %>%
    select(-fam_count)

# Base line: which Ed PGI is most predictive (within families) of Ed years?
feols(edyears ~ 1 + edpgi_all_imputed_self, data = full.data) %>% summary()
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


lm(edpgi_all_imputed_self ~ 1 + edpgi_exclude_imputed_random, data = full.data) %>% summary()



ivreg(as.integer(edyears >= 18) ~ 1 + edpgi_all_imputed_self + edpgi_all_imputed_parental | edpgi_all_imputed_random + edpgi_all_imputed_parental, data = full.data) %>% summary()

ivreg(as.integer(edyears >= 18) ~ 1 + edpgi_all_imputed_self | edpgi_all_imputed_random , data = full.data) %>% summary()


# The gain from ORIV (correlation)
summary(lm(edyears ~ 1 + edpgi_all_imputed_self, data = full.data))
summary(GORIV(edyears ~ 1,
    "edpgi_all_imputed_self", "edpgi_exclude_imputed_self",
    data = full.data, IID = "eid"))

# The gain from ORIV (causal)
source("oriv.R")
summary(lm(edyears ~ 1 + edpgi_all_imputed_self + edpgi_all_imputed_parental,
    data = full.data))
summary(GORIV(edyears ~ 1,
    "edpgi_all_imputed_self", "edpgi_exclude_imputed_self",
    IV.list = c(
        "edpgi_all_imputed_random", "edpgi_exclude_imputed_random"),
    control.list =
        c("edpgi_all_imputed_parental", "edpgi_exclude_imputed_parental"),
    data = full.data, IID = "eid"))

# For higher ed attendance
summary(lm(I(edyears >= 18) ~ 1 + edpgi_all_imputed_self + edpgi_all_imputed_parental,
    data = full.data))
summary(GORIV(I(edyears >= 18) ~ 1,
    "edpgi_all_imputed_self", "edpgi_exclude_imputed_self",
    IV.list = c(
        "edpgi_all_imputed_random", "edpgi_exclude_imputed_random"),
    control.list =
        c("edpgi_all_imputed_parental", "edpgi_exclude_imputed_parental"),
    data = full.data, IID = "eid"))




# Causal, with + with/out ORIV measurement error adjustment.
summary(lm(edyears ~ 1 + edpgi_all_imputed_random, data = full.data))
summary(ivreg(edyears ~ 1 + edpgi_all_imputed_self | edpgi_all_imputed_random , data = full.data))
summary(GORIV(edyears ~ 1,
    "edpgi_all_imputed_random", "edpgi_exclude_imputed_random",
    IID = "eid", FID = "famid",
    data = full.data))
# Causal, with + with/out ORIV measurement error adjustment.
summary(lm(log(soc_median_hourly) ~ 1 + edpgi_all_imputed_random, data = full.data))
summary(GORIV(log(soc_median_hourly) ~ 1,
    "edpgi_all_imputed_random", "edpgi_exclude_imputed_random",
    data = full.data, IID = "eid", FID = "famid"))


################################################################################
## Quasi second-stage, Ed PGI -> Ed Years

# OLS (not causal)
full.data %>%
    lm(edyears ~ 1 + edpgi_self, data = .) %>%
    summary() %>%
    print()

# Show the correlation between the random value, and Ed years
full.data %>%
    lm(edyears ~ 1 + edpgi_random, data = .) %>%
    summary() %>%
    print()

# Show the correlation between the self + parents value, and Ed years
full.data %>%
    lm(edyears ~ 1 + edpgi_self + edpgi_parents, data = .) %>%
    summary() %>%
    print()

# Try it with the random part as an instrument.
full.data %>%
    ivreg::ivreg(edyears ~ 1 + edpgi_self | 1 + edpgi_random,
        data = .) %>%
    summary() %>%
    print()

# Show the mean Ed PGI random component across the parent dist
full.data$parental_quantile <-
    ecdf(full.data$edpgi_parents)(full.data$edpgi_parents)
full.data$parental_quantile <- factor(
    round(10 * full.data$parental_quantile))
# Compare this to Houmark+ (2024) Figure 1.
full.data %>%
    group_by(parental_quantile) %>%
    summarise(
        edpgi_random_mean   = mean(edpgi_random),
        edpgi_random_sd     = sd(edpgi_random),
        edpgi_self_mean     = mean(edpgi_self),
        edpgi_self_sd       = sd(edpgi_self),
        edpgi_parental_mean = mean(edpgi_parents),
        edpgi_parental_sd   = sd(edpgi_parents)) %>%
    View()
mean(full.data$edpgi_self > full.data$edpgi_parents)
# Ensure this does not vary across the parental distribution.
full.data %>%
    lm(edpgi_self ~ 1 + edpgi_random * parental_quantile, data = .) %>%
    summary() %>%
    print()
#  TEST: does this association hold true among people for whom we observe one or both parents? -> Yes.
full.data %>%
    filter(father_present + mother_present > 0) %>%
    lm(edpgi_random ~ 1 + 
        sex_male * (edpgi_father + edpgi_mother), data = .) %>%
    summary() %>%
    print()
# Compare to the Ed PGI in raw form.
full.data %>%
    filter(father_present + mother_present > 0) %>%
    lm(edpgi_self ~ 1 +
        sex_male * (edpgi_father + edpgi_mother), data = .) %>%
    summary() %>%
    print()
