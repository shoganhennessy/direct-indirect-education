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
fig.width <- 7.5
fig.height <- fig.width
presentation.width <- (3 / 2) * fig.width
presentation.height <- (2 / 3) * presentation.width
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
## Plot: correlates with the Ed PGI, and whether they are for random component.

# Show how Ed PGI is related to demographics, lm(Z ~ 1 + X)


# Get the demographic possible confounders.
demographic.data <- analysis.data %>%
    select(
        edpgi_all_imputed_self,
        edpgi_all_imputed_random,
        sex_male,
        recruitedage,
        sibling_count,
        #TODO: urban,
        #TODO: adhd_pgi,
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
    "t2diabetes_pgi"    = "Other PGI: Diabetes (T2)",
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
        name = TeX(r"(Association Estimates $Z_i = \beta\,' \, X_i + \epsilon_i$)"),
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
        #name = TeX(r"(Association Estimates $Z_i = \beta\,' \, X_i + \epsilon_i$)"),
        breaks = seq(-1, 1, by = 0.05)) +
    ggtitle(TeX("Demographic Information")) +
    theme(plot.title = element_text(size = rel(1), hjust = 0),
        plot.title.position = "plot",
        plot.margin = unit(c(0.5, 3, 0, 0), "mm"),
        legend.position = c(0.85, 0.95),
        axis.text.y = element_text(hjust = 0))

# Save this plot
ggsave(file.path(figures.folder, "demographic-correlates.png"),
    plot = combined.plot,
    units = "cm",
    width = 1.25 * presentation.width, height = presentation.height)


################################################################################
## Testing out the first-stage, with random component of Ed PGI (both sources).
# For notation, write Z for main Ed PGI using Okbay (2022) weights (23&Me)
# and Ztilde for the supplementary Ed PGI from UKB subsample.

# Random component Z -> Z
Z_Z_random.reg <- lm(edpgi_all_imputed_self ~ 1 + 
    edpgi_all_imputed_random,# + edpgi_all_imputed_parental,
    data = analysis.data)
Z_Z_random.reg %>% summary() %>% print()

# Random component Ztilde -> Z
Ztilde_Z_random.reg <- lm(edpgi_all_imputed_self ~ 1 + 
    edpgi_exclude_imputed_random,# + edpgi_exclude_imputed_parental,
    data = analysis.data)
Ztilde_Z_random.reg %>% summary() %>% print()

# Random component Z -> Ztilde
Z_Ztilde_random.reg <- lm(edpgi_exclude_imputed_self ~ 1 + 
    edpgi_all_imputed_random,# + edpgi_all_imputed_parental,
    data = analysis.data)
Z_Ztilde_random.reg %>% summary() %>% print()

# Random component Ztilde -> Ztilde
Ztilde_Ztilde_random.reg <- lm(edpgi_exclude_imputed_self ~ 1 +
    edpgi_exclude_imputed_random,# + edpgi_exclude_imputed_parental,
    data = analysis.data)
Ztilde_Ztilde_random.reg %>% summary() %>% print()

## Estimate as an ORIV first-stage (stacked).
source("oriv.R")
# Z main
oriv_Z_random.reg <- GORIV(edpgi_all_imputed_self ~ 1,# +
    #edpgi_all_imputed_parental + edpgi_exclude_imputed_parental,
    "edpgi_all_imputed_random", "edpgi_exclude_imputed_random",
    IID = "eid", FID = "famid",
    data = analysis.data)
oriv_Z_random.reg %>% summary() %>% print()
# Z tilde
oriv_Ztilde_random.reg <- GORIV(edpgi_exclude_imputed_self ~ 1,
    "edpgi_all_imputed_random", "edpgi_exclude_imputed_random",
    IID = "eid", FID = "famid",
    data = analysis.data)
oriv_Ztilde_random.reg %>% summary() %>% print()

#TODO: Adjust the ORIV to have PGI^random as the IV 
#TODO: and the PGI as the treatment (thus 4 components instead of 2).

#TODO: Put everything into a LaTeX table (see notes for shape).
#TODO: Adjust the UKB code so that chr_ data are phased.

################################################################################
## LaTeX Table of the genetic first-stage.

# Compose all the point estimates and SEs, then F stat, obs count, R^2.
point_est.data <- data.frame(
    # column 1: Z -> Z
    col1 = c(
        coef(summary(Z_Z_random.reg))["edpgi_all_imputed_random", "Estimate"],
        coef(summary(Z_Z_random.reg))["edpgi_all_imputed_random", "Std. Error"],
        NA, NA, NA, NA,
        round(as.numeric(summary(Z_Z_random.reg)$fstatistic["value"])),
        summary(Z_Z_random.reg)$r.squared,
        NROW(analysis.data)),
    # column 2: Ztilde -> Z
    col2 = c(NA, NA,
        coef(summary(Ztilde_Z_random.reg))["edpgi_exclude_imputed_random", "Estimate"],
        coef(summary(Ztilde_Z_random.reg))["edpgi_exclude_imputed_random", "Std. Error"],
        NA, NA,
        round(as.numeric(summary(Ztilde_Z_random.reg)$fstatistic["value"])),
        summary(Ztilde_Z_random.reg)$r.squared,
        NROW(analysis.data)),
    # column 3: obviously related Z, Ztilde -> Z
    col3 = c(NA, NA, NA, NA,
        coeftable(oriv_Z_random.reg)["fit_PGI_MAIN", "Estimate"],
        coeftable(oriv_Z_random.reg)["fit_PGI_MAIN", "Std. Error"],
        round(fitstat(oriv_Z_random.reg, "ivf")[1]$`ivf1::PGI_MAIN`$stat),
        NA,
        NROW(analysis.data)),
    # Column 4: Z -> Ztilde
    col4 = c(
        coef(summary(Z_Ztilde_random.reg))["edpgi_all_imputed_random", "Estimate"],
        coef(summary(Z_Ztilde_random.reg))["edpgi_all_imputed_random", "Std. Error"],
        NA, NA, NA, NA,
        round(as.numeric(summary(Z_Ztilde_random.reg)$fstatistic["value"])),
        summary(Z_Ztilde_random.reg)$r.squared,
        NROW(analysis.data)),
    # column 5: Ztilde -> Ztilde
    col5 = c(NA, NA,
        coef(summary(Ztilde_Ztilde_random.reg))["edpgi_exclude_imputed_random", "Estimate"],
        coef(summary(Ztilde_Ztilde_random.reg))["edpgi_exclude_imputed_random", "Std. Error"],
        NA, NA,
        round(as.numeric(summary(Ztilde_Ztilde_random.reg)$fstatistic["value"])),
        summary(Ztilde_Ztilde_random.reg)$r.squared,
        NROW(analysis.data)),
    # column 6: obviously related Z, Ztilde -> Z
    col6 = c(NA, NA, NA, NA,
        coeftable(oriv_Ztilde_random.reg)["fit_PGI_MAIN", "Estimate"],
        coeftable(oriv_Ztilde_random.reg)["fit_PGI_MAIN", "Std. Error"],
        round(fitstat(oriv_Ztilde_random.reg, "ivf")[1]$`ivf1::PGI_MAIN`$stat),
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
point_est.data$names = c("Primary Ed PGI", " ",
        "Secondary Ed PGI", " ",
        "Combined (Obvious Related)", " ",
        "\\\\[-1.8ex]\\hline \\\\[-1.8ex] $F$-statistics",
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
summary(lm(edyears ~ 1 + edpgi_all_imputed_self, data = analysis.data))
summary(GORIV(edyears ~ 1,
    "edpgi_all_imputed_self", "edpgi_exclude_imputed_self",
    data = analysis.data, IID = "eid", FID = "famid"))
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
