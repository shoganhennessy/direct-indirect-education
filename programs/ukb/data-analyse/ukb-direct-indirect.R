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
## 1. Mediation, taking education returns as correlation.

# Define a function, which runs the mediation analysis.
outcome.list <- c("log_soc_median_hourly",
    "log_soc_median_annual", "log_householdincome_midpoint")
regressors.entry <- "~ 1 + edpgi_all_imputed_self + edpgi_all_imputed_parental"
mediator.entry <- "~ 1 + edyears * edpgi_all_imputed_self + edpgi_all_imputed_parental"
control.list <- paste0(
    "+ adhd_pgi + asthma_pgi + bipolar_pgi + bmi_pgi + height_pgi + schizophrenia_pgi + t2diabetes_pgi",
    "+ visityear + sibling_count + sex_male + birthyear + urban + mother_present + father_present")

# Define a function, which runs the mediation analysis.
mediate.table <- function(outcome.entry, control.entry, given.data,
    boot.samples = 10^3) {
    # Get the relevant data.
    reg.data <- given.data[!is.na(given.data[[outcome.entry]]), ]
    # Total effect.
    total.reg <- lm(
        formula(paste0(outcome.entry, regressors.entry, control.entry)),
        data = reg.data)
    print(summary(total.reg))
    # Mediation first-stage
    firststage.reg <- lm(
        formula(paste0("edyears", regressors.entry, control.entry)),
        data = reg.data)
    print(summary(firststage.reg))
    # Mediation second-stage
    secondstage.reg <- lm(
        formula(paste0(outcome.entry, mediator.entry, control.entry)),
        data = reg.data)
    print(summary(secondstage.reg))
    # Get implied education returns.
    ed.returns.reg <- avg_slopes(secondstage.reg, hypothesis = "edyears = 0")
    # Estimate the mediation.
    mediation.reg <- mediate(firststage.reg, secondstage.reg,
        data = analysis.data,
        treat = "edpgi_all_imputed_self", mediator = "edyears",
        robustSE = FALSE, sims = boot.samples)
    print(summary(mediation.reg))
    # Get the effect estimates.
    firststage.point <- summary(firststage.reg)$coefficients["edpgi_all_imputed_self", "Estimate"]
    firststage.se <- summary(firststage.reg)$coefficients["edpgi_all_imputed_self", "Std. Error"]
    total.point <- summary(total.reg)$coefficients["edpgi_all_imputed_self", "Estimate"]
    total.se <- summary(total.reg)$coefficients["edpgi_all_imputed_self", "Std. Error"]
    ed.returns.point <- ed.returns.reg$estimate
    ed.returns.se <- ed.returns.reg$std.error
    # Get the mediation estimates.
    direct.point <- mediation.reg$z.avg
    direct.se <- sd(mediation.reg$z.avg.sims)
    indirect.point <- mediation.reg$d.avg
    indirect.se <- sd(mediation.reg$d.avg.sims)
    percent.point <- mediation.reg$n.avg
    percent.se <- sd(mediation.reg$n.avg.sims)
    obs.count <- mediation.reg$nobs
    # Give a table at the end.
    results <- data.frame(value = c(
        sprintf("%.3f",   firststage.point),
        sprintf("(%.3f)", firststage.se),
        sprintf("%.3f",   total.point),
        sprintf("(%.3f)", total.se),
        sprintf("%.3f",   ed.returns.point),
        sprintf("(%.3f)", ed.returns.se),
        sprintf("%.3f",   direct.point),
        sprintf("(%.3f)", direct.se),
        sprintf("%.3f",   indirect.point),
        sprintf("(%.3f)", indirect.se),
        sprintf("%.3f",   percent.point),
        sprintf("(%.3f)", percent.se),
        prettyNum(obs.count, big.mark = ",", scientific = FALSE)))
    return(results)
}

# Define the names of the table.
table.entries <- data.frame(value = c(
    "Education Effect", "",
    "Total Genetic Effect", "",
    "Education Returns", "",
    r"(\midrule Direct, Ed PGI effect)", "",
    "Indirect, Education Effect", "",
    "Percent mediated through education", "",
    "Observation count"))

# Loop across the outcomes.
for (outcome.entry in outcome.list){
    print(outcome.entry)
    # Estimate on the outcome (with controls).
    controls.est <- mediate.table(outcome.entry, "", analysis.data,
        boot.samples = 10^3)
    # Save the entry
    table.entries <- cbind(table.entries, controls.est)
    # Estimate on the outcome (with controls).
    nocontrols.est <- mediate.table(outcome.entry, control.list, analysis.data,
        boot.samples = 10^3)
    # Save the entry
    table.entries <- cbind(table.entries, nocontrols.est)
}
print(table.entries)

# Save the LaTeX table
rbind(table.entries[1:12, ],
    c(r"(\midrule Controls included?)", rep(c("No", "Yes"), 3)),
    table.entries[13, ]) %>%
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
        file = file.path(tables.folder, "mediate-correlational.tex"))


################################################################################
## Sensitivity analysis, for different education returns values.

mediate.sensitivity <- function(beta, outcome.entry, given.data, indices, boot = FALSE) {
    # Get the relevant data.
    if (boot == TRUE){
        reg.data <- given.data[indices, ]
    } else{
        reg.data <- given.data
    }
    # Total effect.
    total.reg <- lm(
        formula(paste0(outcome.entry, regressors.entry, control.list)),
        data = reg.data)
    # Mediation first-stage
    firststage.reg <- lm(
        formula(paste0("edyears", regressors.entry, control.list)),
        data = reg.data)
    # Collect the mediation est, with given beta value.
    indirect.effect <- (
        coeftable(firststage.reg)["edpgi_all_imputed_self", "Estimate"] * beta)
    total.effect <- coeftable(total.reg)["edpgi_all_imputed_self", "Estimate"]
    direct.effect <- total.effect - indirect.effect
    # Give the entry
    return(c(
        total.effect         = total.effect,
        direct.effect        = direct.effect,
        indirect.effect      = indirect.effect,
        indirect.share       = indirect.effect / total.effect))
}

# Estimate the Ses by the bootstrap.
mediate.sensitivity.boot <- function(beta.entry, outcome.entry, given.data, R = 10^3) {
    # Point estimates from the original data
    point.est <- mediate.sensitivity(beta.entry, outcome.entry, given.data,
        indices = 1:nrow(given.data), boot = FALSE)
    # Boot-compatible wrapper that calls existing function
    boot.wrapper <- function(data, indices) {
        mediate.sensitivity(beta.entry, outcome.entry, data, indices, boot = TRUE)
    }
    # Run the bootstrap
    boot.out <- boot(data = given.data, statistic = boot.wrapper, R = R)
    # Extract SEs as the SD of the bootstrap distribution for each statistic
    boot.ses <- apply(boot.out$t, 2, sd, na.rm = TRUE)
    names(boot.ses) <- paste0(names(point.est), ".se")
    # Return as a one-row dataframe
    result <- as.data.frame(t(c(beta = beta.entry, point.est, boot.ses)))
    return(result)
}
# Loop across relevant values.
beta.list <- seq(0, 0.2, by = 0.025)
library(boot)
sensitivity.data <- do.call(rbind, lapply(beta.list, function(beta.entry) {
    print(beta.entry)
    mediate.sensitivity.boot(beta.entry, outcome.list[1], analysis.data, R = 10^3)
}))
print(sensitivity.data)

# Get the correlational estimate for returns to education.
edreturns.reg <- lm(formula(paste0(
    outcome.list[1], mediator.entry, control.list)), data = analysis.data)
ed.returns.point <- avg_slopes(edreturns.reg, hypothesis = "edyears = 0")$estimate

# Define the sensitivity plot, for direct and indirect channels.
direct.plot <- sensitivity.data %>%
    ggplot(aes(x = 100 * beta)) +
    geom_vline(xintercept = 100 * ed.returns.point, linetype = "dashed", alpha = 0.5) +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.2) +
    annotate("text", x = 4.5, y = 0.1, label = "OLS") +
    # 1. Direct effect
    geom_line(aes(y = direct.effect), colour = colour.list[1], linewidth = 1) +
    geom_ribbon(fill = colour.list[1], alpha = 0.2, aes(
        ymin = direct.effect - 1.96 * direct.effect.se,
        ymax = direct.effect + 1.96 * direct.effect.se)) +
    annotate("text", colour = colour.list[1],
        x = 12.5, y = -0.05, fontface = "bold", label = "Direct") +
    # 2. Indirect effect
    geom_line(aes(y = indirect.effect), colour = colour.list[2], linewidth = 1) +
    geom_ribbon(fill = colour.list[2], alpha = 0.2, aes(
        ymin = indirect.effect - 1.96 * indirect.effect.se,
        ymax = indirect.effect + 1.96 * indirect.effect.se)) +
    annotate("text", colour = colour.list[2],
        x = 12.5, y = 0.1, fontface = "bold", label = "Indirect") +
    # 3. Total effect
    geom_line(aes(y = total.effect), colour = colour.list[3], linewidth = 1) +
    geom_ribbon(fill = colour.list[3], alpha = 0.2, aes(
        ymin = total.effect - 1.96 * total.effect.se,
        ymax = total.effect + 1.96 * total.effect.se)) +
    annotate("text", colour = colour.list[3],
        x = 12.5, y = 0.025, fontface = "bold", label = "Total") +
    theme_bw() +
    scale_x_continuous(expand = c(0, 0),
        name = "Education Return Value, Percent",
        limits = 100.5 * c(min(beta.list), max(beta.list)),
        breaks = 100 * seq(min(beta.list), max(beta.list), by = 0.05)) +
    scale_y_continuous(expand = c(0, 0),
        name = "",
        limits = c(-0.076, 0.13),
        breaks = seq(-0.1, 0.5, by = 0.025)) +
    ggtitle("Direct and Indirect Effect Estimates") +
    theme(plot.title = element_text(hjust = 0, size = rel(1)),
        plot.title.position = "plot",
        plot.margin = unit(c(0.5, 3, 1.5, 0.25), "mm"),
        axis.title.x = element_text(vjust = -0.25))
# Save the plot.
ggsave(file.path(figures.folder, "sens-hourly-wage.png"),
    plot = direct.plot,
    units = "cm", dpi = 300, width = 2 * fig.width, height = fig.height)


################################################################################
## Sensitivity analysis, all three outcomes, as a percent figure.

# SHow the relevant outcome variables
print(outcome.list)

# Operate the sensitivity analysis for the first.
wage_sensitivity.data <- sensitivity.data
# Operate the sensitivity analysis for the second
income_sensitivity.data <- do.call(rbind, lapply(beta.list, function(beta.entry) {
    print(beta.entry)
    mediate.sensitivity.boot(beta.entry, outcome.list[2], analysis.data, R = 10^3)
}))
# Operate the sensitivity analysis for the third
house_sensitivity.data <- do.call(rbind, lapply(beta.list, function(beta.entry) {
    print(beta.entry)
    mediate.sensitivity.boot(beta.entry, outcome.list[3], analysis.data, R = 10^3)
}))
# Combine these data, to plot.
sensitivity_combined.data <- rbind(
    mutate(wage_sensitivity.data, outcome = "Occ hourly wage"),
    mutate(income_sensitivity.data, outcome = "Occ annual income"),
    mutate(house_sensitivity.data, outcome = "Occ household income (midpoint imputed)"))
print(sensitivity_combined.data)

# Define the sensitivity plot, for percent through education
percent.plot <- sensitivity_combined.data %>%
    ggplot(aes(x = 100 * beta)) +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.2) +
    geom_hline(yintercept = 1, linetype = "dashed", alpha = 0.2) +
    # plot percent, by different outcome
    geom_line(aes(y = indirect.share), colour = "orange", linewidth = 1) +
    geom_ribbon(alpha = 0.2, fill = "orange", aes(
        ymin = indirect.share - 1.96 * indirect.share.se,
        ymax = indirect.share + 1.96 * indirect.share.se)) +
    theme_bw() +
    scale_x_continuous(expand = c(0.05, 0),
        name = "Education Return Value, Percent",
        limits = 100.5 * c(min(beta.list), max(beta.list)),
        breaks = 100 * seq(min(beta.list), max(beta.list), by = 0.05)) +
    scale_y_continuous(expand = c(0, 0),
        name = "",
        breaks = seq(-1, 2, by = 0.25)) +
        coord_cartesian(ylim = c(-0.1, 1.55)) + 
    ggtitle("Indirect Effect / Total Effect") +
    theme(plot.title = element_text(hjust = 0, size = rel(1)),
        plot.title.position = "plot",
        plot.margin = unit(c(0.5, 3, 1.5, 0.25), "mm"),
        axis.title.x = element_text(vjust = -0.25)) +
    facet_wrap(~ outcome)
# Save the plot.
ggsave(file.path(figures.folder, "sens-percent.png"),
    plot = percent.plot,
    units = "cm", dpi = 300,
    width = 3 * fig.width, height = fig.height)


################################################################################
## Sensitivity analysis, based on prior literature estimates.

# Get the external IV estimates of education returns.
edreturns.data <- data.folder %>%
    file.path("..", "..", "ed-returns", "causal-returns-database.xlsx") %>%
    readxl::read_excel(sheet = 2)
# Get the British IV estimates.
beta.dist <- edreturns.data %>%
    filter(Country == "UK") %>%
    mutate(causalest = IV / 100) %>%
    pull(causalest)

# Get the indices to sample from.
indices.list <- 1:nrow(analysis.data)

# Bootstrap mediation estimates for each outcome.
boot.sims <- 10^4
wage_boot.data <- read.csv(text =
    "beta,total.effect,direct.effect,indirect.effect,indirect.share")
income_boot.data <- wage_boot.data
house_boot.data <- wage_boot.data

# Loop across bootstrap distribution.
for (i in 1:boot.sims){
    if ((100 * (i / boot.sims)) %% 1 == 0) {
        print(paste0(i, " out of ", boot.sims, ": ", 100 * i / boot.sims, "%."))
    }
    # Calculate a bootstrap entry.
    beta.boot <- sample(beta.dist, 1)
    indices.boot <- sample(indices.list, replace = TRUE)
    # 1. Calculate a wage bootstrap entry
    wage.boot <- mediate.sensitivity(beta.boot, outcome.list[1], analysis.data,
        indices.boot, boot = TRUE)
    wage_boot.data <- rbind(wage_boot.data,
        data.frame(t(c(beta = beta.boot, wage.boot))))
    # 2. Calculate an income bootstrap entry
    income.boot <- mediate.sensitivity(beta.boot, outcome.list[2], analysis.data,
        indices.boot, boot = TRUE)
    income_boot.data <- rbind(income_boot.data,
        data.frame(t(c(beta = beta.boot, income.boot))))
    # 3. Calculate a household income bootstrap entry
    house.boot <- mediate.sensitivity(beta.boot, outcome.list[3], analysis.data,
        indices.boot, boot = TRUE)
    house_boot.data <- rbind(house_boot.data,
        data.frame(t(c(beta = beta.boot, house.boot))))
}


# Define a function to take the boot data, and given a column with.
boot.table <- function(boot.data, outcome.entry) {
    # Get the effect estimates.
    total.point      <- mean(boot.data$total.effect)
    total.se         <- sd(boot.data$total.effect)
    ed.returns.point <- mean(boot.data$beta)
    ed.returns.se    <- sd(boot.data$beta)
    # Get the mediation estimates.
    direct.point   <- mean(boot.data$direct.effect)
    direct.se      <- sd(boot.data$direct.effect)
    indirect.point <- mean(boot.data$indirect.effect)
    indirect.se    <- sd(boot.data$indirect.effect)
    percent.point  <- mean(boot.data$indirect.share)
    percent.se     <- sd(boot.data$indirect.share)
    obs.count      <- nrow(analysis.data[!is.na(analysis.data[[outcome.entry]]), ])
    # Give a table at the end.
    results <- data.frame(value = c(
        sprintf("%.3f",   total.point),
        sprintf("(%.3f)", total.se),
        sprintf("%.3f",   ed.returns.point),
        sprintf("(%.3f)", ed.returns.se),
        sprintf("%.3f",   direct.point),
        sprintf("(%.3f)", direct.se),
        sprintf("%.3f",   indirect.point),
        sprintf("(%.3f)", indirect.se),
        sprintf("%.3f",   percent.point),
        sprintf("(%.3f)", percent.se),
        prettyNum(obs.count, big.mark = ",", scientific = FALSE)))
    return(results)
}

# Define the names of the table.
boot_table.data <- data.frame(value = c(
    "Total Genetic Effect", "",
    "Education Returns", "",
    r"(\midrule Direct, Ed PGI effect)", "",
    "Indirect, Education Effect", "",
    "Percent mediated through education", "",
    r"(\midrule Observation count)")) %>%
    cbind(
        boot.table(wage_boot.data, outcome.list[1]),
        boot.table(income_boot.data, outcome.list[2]),
        boot.table(house_boot.data, outcome.list[3]))

# Save the LaTeX table.
boot_table.data %>%
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
        file = file.path(tables.folder, "mediate-boot.tex"))
