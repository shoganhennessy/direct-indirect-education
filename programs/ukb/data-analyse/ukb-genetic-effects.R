#!/usr/bin/R
## Senan Hogan-Hennessy, 24 March 2025
## UKB data -> Genetic effects, Ed PGI random component -> Ed Years + Income.
set.seed(47)
print(Sys.time())

# Functions for data manipulation and visualisation
library(tidyverse)
library(gtsummary)
# Functions for tables into TeX
library(xtable)
# The standard, linear, IV + FE estimator package.
library(fixest)
# Define number of digits in tables and graphs
digits.no <- 2
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
## Convenience functions.

# Define a function to automate Bin scatter.
binscatter.plot <- function(data, x, y, colour.name,
    option = "none", half.line.slope = 1, half.line.intercept = FALSE){
    library(binsreg)
    # Run the binscatter regression.
    binscatter.data <- binsreg(data[[y]], data[[x]], randcut = 1,
            line = c(1, 1), ci = c(1, 1), cb = c(1, 1), polyreg = 3)  %>%
        magrittr::use_series("data.plot") %>%
        magrittr::use_series("Group Full Sample")
    # Plot the binscatter
    binscatter.ggplot <- ggplot() +
        # Add the bin means
        geom_point(data = binscatter.data$data.dots, aes(x = x, y = fit),
            colour = colour.name, size = 2) +
        # Add the bin error bars
        #geom_errorbar(data = binscatter.data$data.ci,
        #    aes(x = x, ymin = ci.l, ymax = ci.r),
        #    colour = colour.name, size = 0.5, width = 0.02, linetype = "solid") +
        geom_ribbon(data = binscatter.data$data.cb,
            aes(x = x, ymin = cb.l, ymax = cb.r),
            fill = colour.name, alpha = 0.2) +
        # Add the line between bin means
        geom_smooth(data = binscatter.data$data.dots, aes(x = x, y = fit),
            se = FALSE, colour = colour.name, size = 0.5, linetype = "solid")
        if (option == "half-line"){
            print("half-line")
            if (half.line.intercept == FALSE){
                half.line.intercept <- median(binscatter.data$data.dots$fit, na.rm = TRUE)
            }
            binscatter.ggplot <- binscatter.ggplot +
                geom_smooth(data = binscatter.data$data.dots,
                    aes(x = x, y = (half.line.intercept + fit * half.line.slope)),
                        se = FALSE, colour = "orange", size = 1, linetype = "solid")
        }
    # Return the ggplot of this.
    return(binscatter.ggplot)
}

# Load my coded version of ORIV.
#source("oriv.R")


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

# Show the construction
print(analysis.data)
print(names(analysis.data))


################################################################################
## Effect on education years Ed PGI -> education, income.

# Estimate effect of PGI on outcomes.
edyears.reg <- feols(edyears
    ~ 1 + edpgi_all_imputed_self + edpgi_all_imputed_parental
    + adhd_pgi + asthma_pgi + bipolar_pgi + bmi_pgi + height_pgi + schizophrenia_pgi + t2diabetes_pgi
    | visityear + sex_male + birthyear + urban + sibling_count + mother_present + father_present,
    data = analysis.data)
print(summary(edyears.reg))

# Define the same piece-wise
outcome.list <- c("edyears", "agefinishededuc", "edqual_highered",
    "log_soc_median_hourly",
    "log_soc_median_annual", "log_householdincome_midpoint")
regressors.entry <- "~ 1 + edpgi_all_imputed_self + edpgi_all_imputed_parental"
control.list <- paste0(
    "+ adhd_pgi + asthma_pgi + bipolar_pgi + bmi_pgi + height_pgi + schizophrenia_pgi + t2diabetes_pgi",
    " | visityear + sibling_count + sex_male + birthyear + urban + mother_present + father_present")

# Run the regression across a list.
reg.list <- list()
i <- 0
for (outcome.entry in outcome.list){
    print(outcome.entry)
    i <- i + 1
    reg.list[[i]] <- feols(formula(
        paste0(outcome.entry, regressors.entry, control.list)),
        data = analysis.data)
}

# Define a function to get the table of coefficients extracted.
extract.table <- function(model, data, outcome.name, causal = TRUE) {
    # Coefficient table
    table <- summary(model)$coeftable
    # Extract estimates and standard errors
    beta_self  <- table["edpgi_all_imputed_self", "Estimate"]
    se_self    <- table["edpgi_all_imputed_self", "Std. Error"]
    if (causal == TRUE){
        beta_parent <- table["edpgi_all_imputed_parental", "Estimate"]
        se_parent   <- table["edpgi_all_imputed_parental", "Std. Error"]
    } else {
        beta_parent <- NA
        se_parent <- NA
    }
    # Model statistics
    r2  <- summary(model)$sq.cor
    n   <- nobs(model)
    outcome_mean <- mean(data[[outcome.name]], na.rm = TRUE)
    if ("log(" %in% outcome.name){
        outcome_mean <- mean(exp(data[[outcome.name]]), na.rm = TRUE)
    }
    # Construct formatted table
    results <- data.frame(
        value = c(
            sprintf("%.3f", beta_self),
            sprintf("(%.3f)", se_self),
            sprintf("%.3f", beta_parent),
            sprintf("(%.3f)", se_parent),
            sprintf("%.3f", outcome_mean),
            sprintf("%.3f", r2),
            prettyNum(n, big.mark = ",", scientific = FALSE)))
    return(results)
}

# Extract the coefficients
row.term <- c(
    "Ed PGI", "",
    "Parental Ed PGI", "",
    r"( \\[-1.8ex]\hline \\[-1.8ex] Outcome mean)",
    r"(Adjusted $R^2$)",
    "Observations")
results.list <- lapply(seq_along(reg.list), function(i) {
    extract.table(
        model = reg.list[[i]],
        data = analysis.data,
        outcome.name = outcome.list[i])})
results.table <- cbind(row.term, do.call(cbind, results.list))
print(results.table)
# Save the LaTeX table
results.table %>%
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
        file = file.path(tables.folder, "genetic-effects.tex"))


################################################################################
## Estimate the same table with (1) raw OLS (2) Ed PGI no controls (3) fam FEs.

regressors_ols.entry <- "~ 1 + edpgi_all_imputed_self"
no_control.list <- ""
fam_control.list <- paste0(
    "+ adhd_pgi + asthma_pgi + bipolar_pgi + bmi_pgi + height_pgi + schizophrenia_pgi + t2diabetes_pgi",
    " | famid + visityear + sibling_count + sex_male + birthyear + urban + mother_present + father_present")

# Run the regression across a list.
reg_ols.list <- list()
reg_nocontrol.list <- list()
reg_fam.list <- list()
i <- 0
for (outcome.entry in outcome.list){
    print(outcome.entry)
    i <- i + 1
    reg_ols.list[[i]] <- feols(formula(
        paste0(outcome.entry, regressors_ols.entry, no_control.list)),
        data = analysis.data)
    reg_nocontrol.list[[i]] <- feols(formula(
        paste0(outcome.entry, regressors.entry, no_control.list)),
        data = analysis.data)
    reg_fam.list[[i]] <- feols(formula(
        paste0(outcome.entry, regressors_ols.entry, fam_control.list)),
        data = analysis.data)
}

# Define a function to get the table of coefficients extracted.
extract_robust.table <- function(model, data, outcome.name, causal = TRUE) {
    # Coefficient table
    table <- summary(model)$coeftable
    # Extract estimates and standard errors
    beta_self  <- table["edpgi_all_imputed_self", "Estimate"]
    se_self    <- table["edpgi_all_imputed_self", "Std. Error"]
    if (causal == TRUE){
        beta_parent <- table["edpgi_all_imputed_parental", "Estimate"]
        se_parent   <- table["edpgi_all_imputed_parental", "Std. Error"]
    } else{
        beta_parent <- NA
        se_parent <- NA
    }
    # Construct formatted table
    results <- data.frame(
        value = c(
            sprintf("%.3f", beta_self),
            sprintf("(%.3f)", se_self),
            sprintf("%.3f", beta_parent),
            sprintf("(%.3f)", se_parent)))
    return(results)
}

## 1. OLS table.
row_lone.term <- c(
    "Ed PGI", "",
    "Parental Ed PGI", "")
results_ols.list <- lapply(seq_along(reg_ols.list), function(i) {
    extract_robust.table(causal = FALSE,
        model = reg_ols.list[[i]],
        data = analysis.data,
        outcome.name = outcome.list[i])})
results_ols.table <- cbind(row_lone.term, do.call(cbind, results_ols.list)) %>% head(2)
print(results_ols.table)
# Save the LaTeX table
results_ols.table %>%
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
        file = file.path(tables.folder, "genetic-ols.tex"))

# 2. Including parental Ed PGI with no controls.
results_nocontrol.list <- lapply(seq_along(reg_nocontrol.list), function(i) {
    extract_robust.table(causal = TRUE,
        model = reg_nocontrol.list[[i]],
        data = analysis.data,
        outcome.name = outcome.list[i])})
results_nocontrol.table <- cbind(row.term[1:4], do.call(cbind, results_nocontrol.list))
print(results_nocontrol.table)
# Save the LaTeX table
results_nocontrol.table %>%
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
        file = file.path(tables.folder, "genetic-nocontrol.tex"))

# 3. Including parental Ed PGI with fam fixed effects
results_fam.list <- lapply(seq_along(reg_fam.list), function(i) {
    extract.table(causal = FALSE,
        model = reg_fam.list[[i]],
        data = analysis.data,
        outcome.name = outcome.list[i])})
results_fam.table <- cbind(row.term, do.call(cbind, results_fam.list))[-c(3,4, 6), ]
print(results_fam.table)
# Save the LaTeX table
results_fam.table %>%
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
        file = file.path(tables.folder, "genetic-fam.tex"))


################################################################################
## Figure: Ed PGI -> Ed years

# Extract point-estimates from the OLS + ORIV estimates.
edyears.ols <- as.numeric(results_ols.table[1, 2])
edyears.ols.text <- paste(results_ols.table[1, 2], results_ols.table[2, 2])
edyears.causal <- as.numeric(results.table[1, 2])
edyears.causal.text <- paste(results.table[1, 2], results.table[2, 2])

# Show correlation between Ed PGI and edyears in a Bin-scatter plot.
edpgi_edyears.plot <- analysis.data %>%
    binscatter.plot(data = .,
        "edpgi_all_imputed_self", "edyears", colour.list[2],
        option = "half-line",
        half.line.slope = edyears.causal / edyears.ols,
        half.line.intercept = 7) +
    # Annotate OLS
    annotate("text", colour = colour.list[2],
        x = -2.5, y = 17.25,
        fontface = "bold",
        label = paste0("Raw OLS = +", edyears.ols.text),
        size = 4.25, hjust = 0, vjust = 0) +
    # Annotate IV
    annotate("text", colour = "orange",
        x = 0.5, y = 10.25,
        fontface = "bold",
        label = paste0("Causal = +", edyears.causal.text),
        size = 4.25, hjust = 0.5, vjust = 0) +
    annotate("curve", colour = "orange",
        x = 2.125, y = 16,
        xend = 2.375, yend = 15.25,
        linewidth = 0.5,
        curvature = -0.125,
        arrow = arrow(length = unit(0.25, 'cm'))) +
    theme_bw() +
    scale_x_continuous(expand = c(0, 0),
        name = "Ed PGI, s.d. units",
        breaks = seq(-5, 5, by = 1),
        limits = c(-3, 3)) +
    scale_y_continuous(expand = c(0, 0.1),
        name = "",
        limits = c(10, 17.75),
        breaks = seq(0, 20, by = 1)) +
    ggtitle("Education Years") +
    theme(plot.title = element_text(size = rel(1), hjust = 0),
        plot.title.position = "plot",
        plot.margin = unit(c(0.5, 3, 0, 0), "mm"))
# Save this plot
ggsave(file.path(figures.folder, "edpgi-edyears-causal.png"),
    plot = edpgi_edyears.plot,
    units = "cm", width = fig.width, height = fig.height)


################################################################################
## Figure: Ed PGI -> Hourly wages.

# Extract point-estimates from the OLS + causal estimates.
earnings.ols <- as.numeric(results_ols.table[1, 5])
earnings.ols.text <- paste(results_ols.table[1, 5], results_ols.table[2, 5])
earnings.causal <- as.numeric(results.table[1, 5])
earnings.causal.text <- paste(results.table[1, 5], results.table[2, 5])

# Show correlation between Ed PGI and edyears in a Bin-scatter plot.
edpgi_earnings.plot <- analysis.data %>%
    mutate(log_soc_median_hourly = log(soc_median_hourly)) %>%
    binscatter.plot(data = ., "edpgi_all_imputed_self", "soc_median_hourly", colour.list[3],
        option = "half-line",
        half.line.slope = earnings.causal[1] / earnings.ols[1],
        half.line.intercept = 8.75) +
    # Annotate OLS
    annotate("text", colour = colour.list[3],
        x = -2.5, y = 24.25,
        fontface = "bold",
        label = paste0("Raw OLS = +", earnings.ols.text),
        size = 4.25, hjust = 0, vjust = 0) +
    # Annotate IV
    annotate("text", colour = "orange",
        x = 0.5, y = 12.75,
        fontface = "bold",
        label = paste0("Causal = +", earnings.causal.text),
        size = 4.25, hjust = 0.5, vjust = 0) +
    annotate("curve", colour = "orange",
        x = 2, y = 22,
        xend = 2.25, yend = 21,
        linewidth = 0.5,
        curvature = -0.125,
        arrow = arrow(length = unit(0.25, 'cm'))) +
    theme_bw() +
    scale_x_continuous(expand = c(0, 0),
        name = "Ed PGI, s.d. units",
        breaks = seq(-5, 5, by = 1),
        limits = c(-3, 3)) +
    scale_y_continuous(expand = c(0, 0.1),
        name = "",
        limits = c(12, 25.1),
        breaks = seq(0, 50, by = 2.5)) +
    ggtitle("Occupation Hourly Wages, £") +
    theme(plot.title = element_text(size = rel(1), hjust = 0),
        plot.title.position = "plot",
        plot.margin = unit(c(0.5, 3, 0, 0), "mm"))
# Save this plot
ggsave(file.path(figures.folder, "edpgi-earnings-causal.png"),
    plot = edpgi_earnings.plot,
    units = "cm", width = fig.width, height = fig.height)
