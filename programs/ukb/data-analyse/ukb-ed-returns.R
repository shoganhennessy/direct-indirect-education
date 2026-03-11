#!/usr/bin/R
## Senan Hogan-Hennessy, 24 March 2025
## UKB data -> Genetic effects, Ed PGI random component -> Ed Years + Income.
set.seed(47)
print(Sys.time())

# The standard, linear, IV estimator package.
library(fixest)
# Causal medation package, Imai Keele Yamamoto (2010)
library(mediation)
# G-computation of causal effects (allows control interations)
library(marginaleffects)
# Functions for tables into TeX
library(xtable)
# Functions for data manipulation and visualisation
library(tidyverse)
# Define number of digits in tables and graphs
digits.no <- 3
# Size for figures
fig.width <- 7.5
fig.height <- fig.width
presentation.width <- (5 / 3) * fig.width
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
    filter(analysis_sample == 1) %>%
    # calculate log values for convenience
    mutate(
        log_soc_median_hourly = log(soc_median_hourly),
        log_soc_median_annual = log(soc_median_annual),
        log_householdincome_midpoint = log(householdincome_midpoint))


################################################################################
## Returns to education.

# Estimate returns to education.
wages.reg <- feols(log_soc_median_hourly
    ~ 1 + edyears * edpgi_all_imputed_self + edpgi_all_imputed_parental
    + adhd_pgi + asthma_pgi + bipolar_pgi + bmi_pgi + height_pgi + schizophrenia_pgi + t2diabetes_pgi
    | visityear + sex_male + birthyear + urban + sibling_count + mother_present + father_present,
    data = analysis.data)
print(summary(wages.reg))
# Get implied returns to education.
edreturns.coef <- avg_slopes(wages.reg, hypothesis = "edyears = 0")
print(edreturns.coef$estimate)
print(edreturns.coef$std.error)

# Define the same piece-wise
outcome.list <- c("log_soc_median_hourly",
    "log_soc_median_annual", "log_householdincome_midpoint")
regressors.entry <- "~ 1 + edyears"
edpgi.spec <- " * edpgi_all_imputed_self"
parental.spec <- " + edpgi_all_imputed_parental"
control.list <- paste0(
    "+ adhd_pgi + asthma_pgi + bipolar_pgi + bmi_pgi + height_pgi + schizophrenia_pgi + t2diabetes_pgi",
    " | visityear + sibling_count + sex_male + birthyear + urban + mother_present + father_present")

# Run the regression across a list.
ols.list <- list()
control_reg.list <- list()
causal.list <- list()
i <- 0
for (outcome.entry in outcome.list){
    print(outcome.entry)
    i <- i + 1
    ols.list[[i]] <- feols(formula(
        paste0(outcome.entry, regressors.entry, control.list)),
        data = analysis.data)
    control_reg.list[[i]] <- feols(formula(
        paste0(outcome.entry, regressors.entry, edpgi.spec, control.list)),
        data = analysis.data)
    causal.list[[i]] <- feols(formula(
        paste0(outcome.entry, regressors.entry, edpgi.spec,
            parental.spec, control.list)),
        data = analysis.data)
}

# Define a function to get the table of coefficients extracted.
extract.table <- function(model, ed.pgi = TRUE, parental = TRUE) {
    # Coefficient table
    table <- summary(model)$coeftable
    # Extract estimates and standard errors
    beta_self  <- table["edyears", "Estimate"]
    se_self    <- table["edyears", "Std. Error"]
    if (ed.pgi == TRUE){
        beta_edpgi <- table["edpgi_all_imputed_self", "Estimate"]
        se_edpgi   <- table["edpgi_all_imputed_self", "Std. Error"]
        beta_interact <- table["edyears:edpgi_all_imputed_self", "Estimate"]
        se_interact   <- table["edyears:edpgi_all_imputed_self", "Std. Error"]
    } else {
        beta_edpgi <- NA
        se_edpgi <- NA
        beta_interact <- NA
        se_interact <- NA
    }
    if (parental == TRUE){
        beta_parental <- table["edpgi_all_imputed_parental", "Estimate"]
        se_parental   <- table["edpgi_all_imputed_parental", "Std. Error"]
    } else {
        beta_parental <- NA
        se_parental <- NA
    }
    # Model statistics
    edreturns.coef <- avg_slopes(model, hypothesis = "edyears = 0")
    edreturns.est <- edreturns.coef$estimate
    edreturns.se <- edreturns.coef$std.error
    r2  <- summary(model)$sq.cor
    n   <- nobs(model)
    # Construct formatted table
    results <- data.frame(value = c(
        sprintf("%.3f",   beta_self),
        sprintf("(%.3f)", se_self),
        sprintf("%.3f",   beta_edpgi),
        sprintf("(%.3f)", se_edpgi),
        sprintf("%.3f",   beta_interact),
        sprintf("(%.3f)", se_interact),
        sprintf("%.3f",   beta_parental),
        sprintf("(%.3f)", se_parental),
        sprintf("%.3f",   edreturns.est),
        sprintf("(%.3f)", edreturns.se), 
        sprintf("%.3f",   r2),
        prettyNum(n, big.mark = ",", scientific = FALSE)))
    return(results)
}

# Combine all regressions into a table.
item <- extract.table(wages.reg)
cbind(item, item)

row.term <- c(
    "Education years", "",
    "Ed PGI", "",
    r"(Education years $\times$ Ed PGI)", "",
    "Parental Ed PGI", "",
    r"(\midrule Collected Education Returns)", "",
    r"(\midrule $R^2$)",
    "Observation count")
edreturns.table <- cbind(row.term,
    extract.table(ols.list[[1]], ed.pgi = FALSE, parental = FALSE),
    extract.table(control_reg.list[[1]], ed.pgi = TRUE, parental = FALSE),
    extract.table(causal.list[[1]], ed.pgi = TRUE, parental = TRUE),
    extract.table(ols.list[[2]], ed.pgi = FALSE, parental = FALSE),
    extract.table(control_reg.list[[2]], ed.pgi = TRUE, parental = FALSE),
    extract.table(causal.list[[2]], ed.pgi = TRUE, parental = TRUE),
    extract.table(ols.list[[3]], ed.pgi = FALSE, parental = FALSE),
    extract.table(control_reg.list[[3]], ed.pgi = TRUE, parental = FALSE),
    extract.table(causal.list[[3]], ed.pgi = TRUE, parental = TRUE))
edreturns.table[edreturns.table == "(NA)" | edreturns.table == "NA"] <- ""
print(edreturns.table)
# Save the LaTeX table
edreturns.table %>%
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
        file = file.path(tables.folder, "ukb-ed-returns.tex"))
