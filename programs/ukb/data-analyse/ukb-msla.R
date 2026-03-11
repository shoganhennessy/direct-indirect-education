#!/usr/bin/R
## Senan Hogan-Hennessy, 28 Jan 2026.
## UKB data, Effect of mlsa -> Edyears, using modern RDD 
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
library(rdrobust)
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
    "#d62728", # Red
    "Sea Green")      # Dark green (for similar to green variables).
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
full.data <- data.folder %>%
    file.path("ukb-cleaned-pheno.csv") %>%
    read_csv() %>%
    filter(!is.na(edyears) & soc_median_hourly > 0)

# Show the construction
print(full.data)
print(names(full.data))

# Define the birthdate cutoff (born 1957 Sept and later) for mlsa.
cutoff <- 1957.75
# Calculate yearmonth of birth
full.data <- full.data %>%
    mutate(birthyearmonth = birthyear + (birthmonth / 12)) %>%
    mutate(
        mlsa_year = birthyearmonth - cutoff,
        mlsa = as.integer(mlsa_year >= 0),
        mlsa_year = birthyearmonth - cutoff,
        mlsa_year_above = mlsa_year * mlsa,
        mlsa_year_below = mlsa_year * (1 - mlsa),
        soc_median_annual = (
            soc_median_hourly * hours_workweek * 40) / 1000,
        soc_median_annual_all = (
            soc_median_hourly_all * hours_workweek * 40) / 1000,
        # Education outcomes based on qualifications
        gcses_completion = as.integer(edyears >= 11),
        # Education outcomes based on age leaving school.
        agefinishededuc = ifelse(agefinishededuc < 0, NA, agefinishededuc),
        edage16 = as.integer(16 <= agefinishededuc)) %>%
    filter(edyears > 0, soc_median_annual > 0, !is.na(birthyearmonth))

# Show the range of leaving school
full.data %>% pull(agefinishededuc) %>% table(exclude = NULL) %>% print()

# Get the sibling imputed analysis sample.
analysis.data <- full.data %>%
    # Get the sibling imputed analysis sample.
    filter(analysis_sample == 1)


################################################################################
## Estimate the effect of mlsa -> edyears.

## RObust RD Estimates on  UKB data.
# First-stage
summary(rdrobust(x = full.data$birthyearmonth,
    y = full.data$edage16,
    c = cutoff, all = TRUE))
summary(rdrobust(x = full.data$birthyearmonth,
    y = full.data$edyears,
    c = cutoff, all = TRUE))
# Second-stage
summary(rdrobust(x = full.data$birthyearmonth,
    fuzzy = full.data$edage16,
    y = log(full.data$soc_median_hourly),
    c = cutoff, all = TRUE))
summary(rdrobust(x = full.data$birthyearmonth,
    fuzzy = full.data$edyears,
    y = log(full.data$soc_median_hourly),
    c = cutoff, all = TRUE))

#! Test -> age finished school
summary(rdrobust(x = full.data$birthyearmonth,
    fuzzy = full.data$edage16,
    y = log(full.data$soc_median_hourly),
    c = cutoff, all = TRUE))
table(full.data$edyears)
summary(rdrobust(x = full.data$birthyearmonth,
    fuzzy = full.data$edyears,
    y = log(full.data$soc_median_hourly),
    c = cutoff, all = TRUE))

# Compare to OLS.
ols.reg <- feols(log(soc_median_hourly) ~ edyears | birthyear,
    data = full.data)
print(summary(ols.reg))
sibling_fe.reg <- feols(log(soc_median_hourly) ~ edyears | birthyear + famid,
    data = filter(full.data, analysis_sample == 1))
print(summary(sibling_fe.reg))


################################################################################
## Plot the first-stage and the reduced form, with edyears.

# Calculate birth month means.
plot.data <- full.data %>%
    filter(abs(birthyearmonth - cutoff) < 11) %>%
    group_by(birthyearmonth) %>%
    summarise(
        # Education outcomes based on qualifications
        edyears = mean(edyears, na.rm = TRUE),
        gcses_completion = mean(gcses_completion, na.rm = TRUE),
        # Education outcomes based on age leaving school.
        agefinishededuc = mean(agefinishededuc, na.rm = TRUE),
        edage16 = mean(edage16, na.rm = TRUE),
        # Occupation coded income outcomes
        soc_median_hourly = mean(log(soc_median_hourly), na.rm = TRUE),
        soc_median_annual = mean(log(soc_median_annual), na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(
        # Edyears
        edyears_below = ifelse(birthyearmonth < cutoff, edyears, NA),
        edyears_above = ifelse(cutoff <= birthyearmonth, edyears, NA),
        # GCSE completion
        gcses_below = ifelse(birthyearmonth < cutoff, gcses_completion, NA),
        gcses_above = ifelse(cutoff <= birthyearmonth, gcses_completion, NA),
        # Age leaving school
        agefinishededuc_below = ifelse(birthyearmonth < cutoff, agefinishededuc, NA),
        agefinishededuc_above = ifelse(cutoff <= birthyearmonth, agefinishededuc, NA),
        # Finish above age 16
        edage16_below = ifelse(birthyearmonth < cutoff, edage16, NA),
        edage16_above = ifelse(cutoff <= birthyearmonth, edage16, NA),
        # Median wages
        wages_below = ifelse(birthyearmonth < cutoff, soc_median_hourly, NA),
        wages_above = ifelse(cutoff <= birthyearmonth, soc_median_hourly, NA),
        # Median income
        income_below = ifelse(birthyearmonth < cutoff, soc_median_annual, NA),
        income_above = ifelse(cutoff <= birthyearmonth, soc_median_annual, NA))

# Plot predicted edyears, before, after.
edyears.plot <- plot.data %>%
    ggplot(aes(x = birthyearmonth)) +
    # Raw Ed years
    geom_point(aes(y = edyears), colour = colour.list[2], alpha = 0.25) +
    # RDD Ed years
    stat_smooth(aes(y = edyears_below), method = "loess", se = FALSE, colour = colour.list[2], linewidth = 2) +
    stat_smooth(aes(y = edyears_above), method = "loess", se = FALSE, colour = colour.list[2], linewidth = 2) +
    stat_smooth(method = "loess", se = FALSE, aes(y = edyears_above), colour = colour.list[2], linewidth = 2) +
    geom_vline(xintercept = cutoff, linetype = "dashed") +
    theme_bw() +
    scale_x_continuous(expand = c(0, 0),
        name = "Birth Year",
        breaks = seq(1940, 1970, by = 2)) +
    scale_y_continuous(expand = c(0, 0),
        name = "",
        limits = c(12.9, 18.1),
        breaks = seq(0, 18, by = 0.5)) +
    ggtitle(TeX(r"(Mean Education Years)")) +
    theme(plot.title = element_text(size = rel(1), hjust = 0),
        plot.title.position = "plot",
        plot.margin = unit(c(0.5, 3, 1, 0), "mm"))
# Save this plot
ggsave(file.path(figures.folder, "mlsa-edyears-nl.png"),
    plot = edyears.plot,
    units = "cm", width = fig.width, height = fig.height)

# Plot predicted GCSE completion rate , before, after.
gcses.plot <- plot.data %>%
    ggplot(aes(x = birthyearmonth)) +
    # Raw Ed years
    geom_point(aes(y = gcses_completion), colour = colour.list[2], alpha = 0.25) +
    # RDD Ed years
    stat_smooth(aes(y = gcses_below), method = "loess", se = FALSE, colour = colour.list[2], linewidth = 2) +
    stat_smooth(aes(y = gcses_above), method = "loess", se = FALSE, colour = colour.list[2], linewidth = 2) +
    stat_smooth(method = "loess", se = FALSE, aes(y = gcses_above), colour = colour.list[2], linewidth = 2) +
    geom_vline(xintercept = cutoff, linetype = "dashed") +
    theme_bw() +
    scale_x_continuous(expand = c(0, 0),
        name = "Birth Year",
        breaks = seq(1940, 1970, by = 2)) +
    scale_y_continuous(expand = c(0, 0),
        name = "",
        limits = c(0.56, 1.01),
        breaks = seq(0, 1, by = 0.05)) +
    ggtitle(TeX(r"(GCSE completion rate)")) +
    theme(plot.title = element_text(size = rel(1), hjust = 0),
        plot.title.position = "plot",
        plot.margin = unit(c(0.5, 3, 1, 0), "mm"))
# Save this plot
ggsave(file.path(figures.folder, "mlsa-gcses-nl.png"),
    plot = gcses.plot,
    units = "cm", width = fig.width, height = fig.height)

# Plot predicted Finish above age 16 completion, before, after.
edage16.plot <- plot.data %>%
    ggplot(aes(x = birthyearmonth)) +
    # Raw Ed years
    geom_point(aes(y = edage16), colour = colour.list[4], alpha = 0.25) +
    # RDD Ed years
    stat_smooth(method = "loess", se = FALSE, aes(y = edage16_below), colour = colour.list[4], linewidth = 2) +
    stat_smooth(method = "loess", se = FALSE, aes(y = edage16_above), colour = colour.list[4], linewidth = 2) +
    geom_vline(xintercept = cutoff, linetype = "dashed") +
    theme_bw() +
    scale_x_continuous(expand = c(0, 0),
        name = "Birth Year",
        breaks = seq(1940, 1970, by = 2)) +
    scale_y_continuous(expand = c(0, 0),
        name = "",
        limits = c(0.56, 1.01),
        breaks = seq(0, 1, by = 0.05)) +
    ggtitle(TeX(r"(Rate of leaving school aged 16 or older)")) +
    theme(plot.title = element_text(size = rel(1), hjust = 0),
        plot.title.position = "plot",
        plot.margin = unit(c(0, 3, 1, 0), "mm"))
# Save this plot
ggsave(file.path(figures.folder, "mlsa-edage16-nl.png"),
    plot = edage16.plot,
    units = "cm", width = fig.width, height = fig.height)

# Plot predicted age leave school, before, after.
agefinishededuc.plot <- plot.data %>%
    ggplot(aes(x = birthyearmonth)) +
    # Raw Ed years
    geom_point(aes(y = agefinishededuc), colour = colour.list[4], alpha = 0.25) +
    # RDD Ed years
    stat_smooth(method = "loess", se = FALSE, aes(y = agefinishededuc_below), colour = colour.list[4], linewidth = 2) +
    stat_smooth(method = "loess", se = FALSE, aes(y = agefinishededuc_above), colour = colour.list[4], linewidth = 2) +
    geom_vline(xintercept = cutoff, linetype = "dashed") +
    theme_bw() +
    scale_x_continuous(expand = c(0, 0),
        name = "Birth Year",
        breaks = seq(1940, 1970, by = 2)) +
    scale_y_continuous(expand = c(0, 0),
        name = "",
        limits = c(12.9, 18.1),
        breaks = seq(0, 18, by = 0.5)) +
    ggtitle(TeX(r"(Mean age leaving school)")) +
    theme(plot.title = element_text(size = rel(1), hjust = 0),
        plot.title.position = "plot",
        plot.margin = unit(c(0, 3, 1, 0), "mm"))
# Save this plot
ggsave(file.path(figures.folder, "mlsa-agefinishededuc-nl.png"),
    plot = agefinishededuc.plot,
    units = "cm", width = fig.width, height = fig.height)

# Plot wages before and after.
wages.plot <- plot.data %>%
    ggplot(aes(x = birthyearmonth)) +
    # Raw wages
    geom_point(aes(y = exp(soc_median_hourly)), colour = colour.list[3], alpha = 0.25) +
    # RDD wages
    stat_smooth(method = "loess", se = FALSE, aes(y = exp(wages_below)), colour = colour.list[3], linewidth = 2) +
    stat_smooth(method = "loess", se = FALSE, aes(y = exp(wages_above)), colour = colour.list[3], linewidth = 2) +
    geom_vline(xintercept = cutoff, linetype = "dashed") +
    theme_bw() +
    scale_x_continuous(expand = c(0, 0),
        name = "Birth Year",
        breaks = seq(1940, 1970, by = 2)) +
    scale_y_continuous(expand = c(0, 0),
        name = "",
        limits = c(16, 20),
        breaks = seq(0, 30, by = 0.5)) +
    ggtitle("Mean Occupational Hourly Wage, £") +
    theme(plot.title = element_text(size = rel(1), hjust = 0),
        plot.title.position = "plot",
        plot.margin = unit(c(0, 3, 1, 0), "mm"))
# Save this plot
ggsave(file.path(figures.folder, "mlsa-wages-nl.png"),
    plot = wages.plot,
    units = "cm", width = fig.width, height = fig.height)

# Plot annual income before and after.
income.plot <- plot.data %>%
    ggplot(aes(x = birthyearmonth)) +
    # Raw income
    geom_point(aes(y = exp(soc_median_annual)), colour = colour.list[3], alpha = 0.25) +
    # RDD income
    stat_smooth(method = "loess", se = FALSE, aes(y = exp(income_below)), colour = colour.list[3], linewidth = 2) +
    stat_smooth(method = "loess", se = FALSE, aes(y = exp(income_above)), colour = colour.list[3], linewidth = 2) +
    geom_vline(xintercept = cutoff, linetype = "dashed") +
    theme_bw() +
    scale_x_continuous(expand = c(0, 0),
        name = "Birth Year",
        breaks = seq(1940, 1970, by = 2)) +
    scale_y_continuous(expand = c(0, 0),
        name = "",
        limits = c(15.5, 30),
        breaks = seq(0, 30, by = 2)) +
    ggtitle("Mean Occupational Annual Income, £ thousands") +
    theme(plot.title = element_text(size = rel(1), hjust = 0),
        plot.title.position = "plot",
        plot.margin = unit(c(0, 3, 1, 0), "mm"))
# Save this plot
ggsave(file.path(figures.folder, "mlsa-income-nl.png"),
    plot = income.plot,
    units = "cm", width = fig.width, height = fig.height)


################################################################################
## Table of the MSLA first-stage and the reduced form.

## Estimate every model necessary
## Linear -> Qualification deinition
## Poly -> leaving school definition.

# 1. Finish above age 16
linear_edage16.reg <- rdrobust(x = full.data$birthyearmonth,
    y = full.data$gcses_completion, c = cutoff)
print(summary(linear_edage16.reg))
poly_edage16.reg <- rdrobust(x = full.data$birthyearmonth,
    y = full.data$edage16, c = cutoff)
print(summary(poly_edage16.reg))
# 2. Ed years
linear_edyears.reg <- rdrobust(x = full.data$birthyearmonth,
    y = full.data$edyears, c = cutoff)
print(summary(linear_edyears.reg))
poly_edyears.reg <- rdrobust(x = full.data$birthyearmonth,
    y = full.data$agefinishededuc, c = cutoff)
print(summary(poly_edyears.reg))
# 3. Log wages
linear_wages.reg <- rdrobust(x = full.data$birthyearmonth,
    y = log(full.data$soc_median_hourly), c = cutoff)
print(summary(linear_wages.reg))
poly_wages.reg <- linear_wages.reg
print(summary(poly_wages.reg))
# 4. Log income
linear_income.reg <- rdrobust(x = full.data$birthyearmonth,
    y = log(full.data$soc_median_annual), c = cutoff)
print(summary(linear_income.reg))
poly_income.reg <- linear_income.reg
print(summary(poly_income.reg))

# Collect every model into an ordered list.
model.list <- list(
    linear_edage16.reg, poly_edage16.reg,
    linear_edyears.reg, poly_edyears.reg,
    linear_wages.reg,   poly_wages.reg,
    linear_income.reg,  poly_income.reg)
# Get point estimates.
point.est <-   sapply(model.list, function(x) x$coef[3])
point.se <-    sapply(model.list, function(x) x$se[3])
point.fstat <- sapply(model.list, function(x) x$z[3]^2)
point.bw  <-   sapply(model.list, function(x) x$bws[1, 1])
point.obs <-   sapply(model.list, function(x) sum(x$N_h))

# Write to a LaTeX table.
mlsa.table <- style_sigfig(point.est, digits = digits.no, decimals = digits.no) %>%
    rbind(style_sigfig(point.se, digits = digits.no, decimals = digits.no))
mlsa.table[2, ] <- mlsa.table[2, ] %>% paste0("(", ., ")")
mlsa.table <- mlsa.table %>%
    rbind(style_sigfig(point.fstat, digits = digits.no, decimals = digits.no)) %>%
    rbind(style_sigfig(point.bw, digits = digits.no, decimals = digits.no)) %>%
    rbind(style_sigfig(point.obs))
# Add labels for the linear/poly
#mlsa.table <- mlsa.table %>% rbind(rep(c("Yes" , " "), 4))
#mlsa.table <- mlsa.table %>% rbind(rep(c(" ", "Yes"), 4))
# And add obs counts at bottom.
#mlsa.table <- mlsa.table %>% rbind(point.obs)
# Give it names.
mlsa.table <- c("MSLA Effect", " ") %>%
    c(r"(\\ $F$-statistics)",
        r"(CCT bandwidth)", r"(Bandwidth observations)") %>%
    cbind(mlsa.table)

# Show the LaTeX table
mlsa.table %>%
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
        file = file.path(tables.folder, "mlsa-firststage-nl.tex"))


################################################################################
## Returns to Education Table, MSLA second-stage.

## Log Wage: Estimate every model necessary
library(splines)
bandwidth <- 11
tri.dist <- abs((full.data$birthyearmonth - cutoff) / bandwidth)
full.data$kernel.wt <- (1 - tri.dist) * (tri.dist <= 1) / bandwidth

# 1. Finish above age 16, OLS
ols_linear_edage16.reg <- full.data %>%
    feols(log(soc_median_hourly) ~ 1 + mlsa +
        bs(mlsa_year_below) + bs(mlsa_year_above) + gcses_completion,
        weights = full.data$kernel.wt,
        se = "twoway", cluster = c("birthyear", "birthmonth"),
        data = .)
print(summary(ols_linear_edage16.reg))
ols_poly_edage16.reg <- full.data %>%
    feols(log(soc_median_hourly) ~ 1 + mlsa +
        bs(mlsa_year_below) + bs(mlsa_year_above) + edage16,
        weights = full.data$kernel.wt,
        se = "twoway", cluster = c("birthyear", "birthmonth"),
        data = .)
print(summary(ols_poly_edage16.reg))
# 1.5 Finish above age 16, IV
iv_linear_edage16.reg <- rdrobust(x = full.data$birthyearmonth,
    fuzzy = full.data$gcses_completion,
    y = log(full.data$soc_median_hourly), c = cutoff)
print(summary(iv_linear_edage16.reg))
iv_poly_edage16.reg <- rdrobust(x = full.data$birthyearmonth,
    fuzzy = full.data$edage16,
    y = log(full.data$soc_median_hourly), c = cutoff)
print(summary(iv_poly_edage16.reg))

# 2. Ed years
ols_linear_edyears.reg <- full.data %>%
    feols(log(soc_median_hourly) ~ 1 + edyears +
        bs(mlsa_year_below) + bs(mlsa_year_above),
        weights = full.data$kernel.wt,
        se = "twoway", cluster = c("birthyear", "birthmonth"),
        data = .)
print(summary(ols_linear_edyears.reg))
ols_poly_edyears.reg <- full.data %>%
    feols(log(soc_median_hourly) ~ 1 + agefinishededuc +
        bs(mlsa_year_below) + bs(mlsa_year_above),
        weights = full.data$kernel.wt,
        se = "twoway", cluster = c("birthyear", "birthmonth"),
        data = .)
print(summary(ols_poly_edyears.reg))
# 2.5 Finish above age 16, IV
iv_linear_edyears.reg <- rdrobust(x = full.data$birthyearmonth,
    fuzzy = full.data$edyears,
    y = log(full.data$soc_median_hourly), c = cutoff)
print(summary(iv_linear_edyears.reg))
iv_poly_edyears.reg <- rdrobust(x = full.data$birthyearmonth,
    fuzzy = full.data$agefinishededuc,
    y = log(full.data$soc_median_hourly), c = cutoff)
print(summary(iv_poly_edyears.reg))

# Collect every model into an ordered list.
edage16.list <- list(
    ols_linear_edage16.reg, ols_poly_edage16.reg,
    iv_linear_edage16.reg,  iv_poly_edage16.reg)
edyears.list <- list(
    ols_linear_edyears.reg, ols_poly_edyears.reg,
    iv_linear_edyears.reg,  iv_poly_edyears.reg)

# Get point estimates.
edage16.point <- c(
    coeftable(ols_linear_edage16.reg)["gcses_completion", "Estimate"],
    coeftable(ols_poly_edage16.reg)["edage16", "Estimate"],
    iv_linear_edage16.reg$coef[3],
    iv_poly_edage16.reg$coef[3])
edage16.point <- style_sigfig(edage16.point, digits = digits.no, decimals = digits.no)
edage16.point <- c(edage16.point[1:2], rep(" ", 2), edage16.point[3:4], rep(" ", 2))

edyears.point <- c(
    coeftable(ols_linear_edyears.reg)["edyears", "Estimate"],
    coeftable(ols_poly_edyears.reg)["agefinishededuc", "Estimate"],
    iv_linear_edyears.reg$coef[3],
    iv_poly_edyears.reg$coef[3])
edyears.point <- style_sigfig(edyears.point, digits = digits.no, decimals = digits.no)
edyears.point <- c(rep(" ", 2), edyears.point[1:2], rep(" ", 2), edyears.point[3:4])

# Get observation number.
point.obs <- sapply(
    list(ols_linear_edage16.reg, ols_poly_edage16.reg,
        ols_linear_edyears.reg, ols_poly_edyears.reg),
    function(x) format(x$nobs, big.mark = ",", scientific = FALSE)) %>%
    c(sapply(list(iv_linear_edage16.reg, iv_poly_edage16.reg,
        iv_linear_edyears.reg, iv_poly_edyears.reg),
    function(x) format(sum(x$N_h), big.mark = ",", scientific = FALSE)))
# Get SEs
edage16.se <- c(
    coeftable(ols_linear_edage16.reg)["gcses_completion", "Std. Error"],
    coeftable(ols_poly_edage16.reg)["edage16", "Std. Error"],
    iv_linear_edage16.reg$se[3],
    iv_poly_edage16.reg$se[3])
edage16.se <- style_sigfig(edage16.se, digits = digits.no, decimals = digits.no) %>% paste0("(", ., ")")
edage16.se <- c(edage16.se[1:2], rep(" ", 2), edage16.se[3:4], rep(" ", 2))

edyears.se <- c(
    coeftable(ols_linear_edyears.reg)["edyears", "Std. Error"],
    coeftable(ols_poly_edyears.reg)["agefinishededuc", "Std. Error"],
    iv_linear_edyears.reg$se[3],
    iv_poly_edyears.reg$se[3])
edyears.se <- style_sigfig(edyears.se, digits = digits.no, decimals = digits.no) %>% paste0("(", ., ")")
edyears.se <- c(rep(" ", 2), edyears.se[1:2], rep(" ", 2), edyears.se[3:4])

# Get the bandwidth values
table.bw <- c(10, 10,
    10, 10,
    iv_linear_edage16.reg$bws[1, 1],
    iv_poly_edage16.reg$bws[1, 1],
    iv_linear_edyears.reg$bws[1, 1],
    iv_poly_edyears.reg$bws[1, 1]) %>%
    style_sigfig(digits = digits.no, decimals = digits.no)

# Get the observation count
table.nobs <- c(
    ols_linear_edage16.reg$nobs,
    ols_poly_edage16.reg$nobs,
    ols_linear_edyears.reg$nobs,
    ols_poly_edyears.reg$nobs,
    sum(iv_linear_edage16.reg$N_h),
    sum(iv_poly_edage16.reg$N_h),
    sum(iv_linear_edyears.reg$N_h),
    sum(iv_poly_edyears.reg$N_h)) %>%
    format(big.mark = ",", scientific = FALSE)
# Build a LaTeX table.
returns.table <- rbind(edage16.point, edage16.se, edyears.point, edyears.se,
    table.bw, table.nobs)
# Give it names.
returns.table <- c("GCSE completion", " ",
        "Education years", " ",
        r"(\\ CCT bandwidth)", r"(Bandwidth observations)") %>%
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
        file = file.path(tables.folder, "mlsa-secondstage-wage-nl.tex"))

## Log Income: Estimate every model necessary
# 1. Finish above age 16, OLS
ols_linear_edage16.reg <- full.data %>%
    feols(log(soc_median_annual) ~ 1 + mlsa +
        bs(mlsa_year_below) + bs(mlsa_year_above) + gcses_completion,
        weights = full.data$kernel.wt,
        se = "twoway", cluster = c("birthyear", "birthmonth"),
        data = .)
print(summary(ols_linear_edage16.reg))
ols_poly_edage16.reg <- full.data %>%
    feols(log(soc_median_annual) ~ 1 + mlsa +
        bs(mlsa_year_below) + bs(mlsa_year_above) + edage16,
        weights = full.data$kernel.wt,
        se = "twoway", cluster = c("birthyear", "birthmonth"),
        data = .)
print(summary(ols_poly_edage16.reg))
# 1.5 Finish above age 16, IV
iv_linear_edage16.reg <- rdrobust(x = full.data$birthyearmonth,
    fuzzy = full.data$gcses_completion,
    y = log(full.data$soc_median_annual), c = cutoff)
print(summary(iv_linear_edage16.reg))
iv_poly_edage16.reg <- rdrobust(x = full.data$birthyearmonth,
    fuzzy = full.data$edage16,
    y = log(full.data$soc_median_annual), c = cutoff)
print(summary(iv_poly_edage16.reg))

# 2. Ed years
ols_linear_edyears.reg <- full.data %>%
    feols(log(soc_median_annual) ~ 1 + edyears +
        bs(mlsa_year_below) + bs(mlsa_year_above),
        weights = full.data$kernel.wt,
        se = "twoway", cluster = c("birthyear", "birthmonth"),
        data = .)
print(summary(ols_linear_edyears.reg))
ols_poly_edyears.reg <- full.data %>%
    feols(log(soc_median_annual) ~ 1 + agefinishededuc +
        bs(mlsa_year_below) + bs(mlsa_year_above),
        weights = full.data$kernel.wt,
        se = "twoway", cluster = c("birthyear", "birthmonth"),
        data = .)
print(summary(ols_poly_edyears.reg))
# 2.5 Finish above age 16, IV
iv_linear_edyears.reg <- rdrobust(x = full.data$birthyearmonth,
    fuzzy = full.data$edyears,
    y = log(full.data$soc_median_annual), c = cutoff)
print(summary(iv_linear_edyears.reg))
iv_poly_edyears.reg <- rdrobust(x = full.data$birthyearmonth,
    fuzzy = full.data$agefinishededuc,
    y = log(full.data$soc_median_annual), c = cutoff)
print(summary(iv_poly_edyears.reg))

# Collect every model into an ordered list.
edage16.list <- list(
    ols_linear_edage16.reg, ols_poly_edage16.reg,
    iv_linear_edage16.reg,  iv_poly_edage16.reg)
edyears.list <- list(
    ols_linear_edyears.reg, ols_poly_edyears.reg,
    iv_linear_edyears.reg,  iv_poly_edyears.reg)

# Get point estimates.
edage16.point <- c(
    coeftable(ols_linear_edage16.reg)["gcses_completion", "Estimate"],
    coeftable(ols_poly_edage16.reg)["edage16", "Estimate"],
    iv_linear_edage16.reg$coef[3],
    iv_poly_edage16.reg$coef[3])
edage16.point <- style_sigfig(edage16.point, digits = digits.no, decimals = digits.no)
edage16.point <- c(edage16.point[1:2], rep(" ", 2), edage16.point[3:4], rep(" ", 2))

edyears.point <- c(
    coeftable(ols_linear_edyears.reg)["edyears", "Estimate"],
    coeftable(ols_poly_edyears.reg)["agefinishededuc", "Estimate"],
    iv_linear_edyears.reg$coef[3],
    iv_poly_edyears.reg$coef[3])
edyears.point <- style_sigfig(edyears.point, digits = digits.no, decimals = digits.no)
edyears.point <- c(rep(" ", 2), edyears.point[1:2], rep(" ", 2), edyears.point[3:4])

# Get observation number.
point.obs <- sapply(
    list(ols_linear_edage16.reg, ols_poly_edage16.reg,
        ols_linear_edyears.reg, ols_poly_edyears.reg),
    function(x) format(x$nobs, big.mark = ",", scientific = FALSE)) %>%
    c(sapply(list(iv_linear_edage16.reg, iv_poly_edage16.reg,
        iv_linear_edyears.reg, iv_poly_edyears.reg),
    function(x) format(sum(x$N_h), big.mark = ",", scientific = FALSE)))
# Get SEs
edage16.se <- c(
    coeftable(ols_linear_edage16.reg)["gcses_completion", "Std. Error"],
    coeftable(ols_poly_edage16.reg)["edage16", "Std. Error"],
    iv_linear_edage16.reg$se[3],
    iv_poly_edage16.reg$se[3])
edage16.se <- style_sigfig(edage16.se, digits = digits.no, decimals = digits.no) %>% paste0("(", ., ")")
edage16.se <- c(edage16.se[1:2], rep(" ", 2), edage16.se[3:4], rep(" ", 2))

edyears.se <- c(
    coeftable(ols_linear_edyears.reg)["edyears", "Std. Error"],
    coeftable(ols_poly_edyears.reg)["agefinishededuc", "Std. Error"],
    iv_linear_edyears.reg$se[3],
    iv_poly_edyears.reg$se[3])
edyears.se <- style_sigfig(edyears.se, digits = digits.no, decimals = digits.no) %>% paste0("(", ., ")")
edyears.se <- c(rep(" ", 2), edyears.se[1:2], rep(" ", 2), edyears.se[3:4])

# Get the bandwidth values
table.bw <- c(10, 10,
    10, 10,
    iv_linear_edage16.reg$bws[1, 1],
    iv_poly_edage16.reg$bws[1, 1],
    iv_linear_edyears.reg$bws[1, 1],
    iv_poly_edyears.reg$bws[1, 1]) %>%
    style_sigfig(digits = digits.no, decimals = digits.no)

# Get the observation count
table.nobs <- c(
    ols_linear_edage16.reg$nobs,
    ols_poly_edage16.reg$nobs,
    ols_linear_edyears.reg$nobs,
    ols_poly_edyears.reg$nobs,
    sum(iv_linear_edage16.reg$N_h),
    sum(iv_poly_edage16.reg$N_h),
    sum(iv_linear_edyears.reg$N_h),
    sum(iv_poly_edyears.reg$N_h)) %>%
    format(big.mark = ",", scientific = FALSE)
# Build a LaTeX table.
returns.table <- rbind(edage16.point, edage16.se, edyears.point, edyears.se,
    table.bw, table.nobs)
# Add labels for the education qual versus school age leaving.
returns.table <- returns.table %>% rbind(rep(c("Yes" , " "), 4))
returns.table <- returns.table %>% rbind(rep(c(" ", "Yes"), 4))
# Give it names.
returns.table <- c("GCSE completion", " ",
        "Education years", " ",
        r"(\\ CCT bandwidth)", r"(Bandwidth observations)",
        r"(\\[-1.8ex]\hline \\[-1.8ex] Qualification definition)",
        "School age definition") %>%
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
        file = file.path(tables.folder, "mlsa-secondstage-income-nl.tex"))


################################################################################
## Calculate the CCT bandwidth option plot, ed years

# Loop across the bandwidth values.
by.value <- 0.25
max.value <- 10
bw_grid <- seq(2 * by.value, max.value - 2 * by.value, by.value)
ests <- se <- numeric(length(bw_grid))
for (i in seq_along(bw_grid)) {
    print(i)
    fit <- rdrobust(
        x = full.data$birthyearmonth,
        fuzzy = full.data$edyears,
        y = log(full.data$soc_median_hourly),
        c = cutoff,
        h = bw_grid[i])
    ests[i] <- fit$coef[3]
    se[i]   <- fit$se[3]
}

# Get the CCT bandiwdth.
fit <- rdrobust(
    x = full.data$birthyearmonth,
    fuzzy = full.data$edyears,
    y = log(full.data$soc_median_hourly),
    c = cutoff)
cct.bw <- fit$bws[1,1]

# Save the output.
df_plot <- data.frame(bw = bw_grid,
    est = ests,
    lo = ests - 1.96*se,
    hi = ests + 1.96*se)
# Plot the outcome.
bandwidth.plot <- df_plot %>%
    ggplot(aes(x = bw, y = est)) +
    geom_vline(xintercept = cct.bw, linetype = "dotted") +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.25) +
    geom_point(colour = colour.list[2]) +
    geom_line(colour = colour.list[2]) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.2, fill = colour.list[2]) +
    annotate("text", x = cct.bw  * 1.05, y = 0.36,
        label = "CCT \nbandwidth",
        size = 4.25,  hjust = 0, vjust = 1) +
    theme_bw() +
    scale_x_continuous(expand = c(0, 0),
        name = "RD bandwidth choice, years to Sept 1957 birth date",
        limits = c(0, max.value - 2 * by.value),
        breaks = seq(0, max.value, by = 4 * by.value)) +
    scale_y_continuous(expand = c(0, 0),
        name = "",
        limits = c(-0.11, 0.41),
        breaks = seq(-1, 3, by = 0.1)) +
    ggtitle(TeX(r"(Fuzzy RD, Education years on Log hourly wage)")) +
    theme(plot.title = element_text(size = rel(1), hjust = 0),
        plot.title.position = "plot",
        plot.margin = unit(c(0.5, 3, 1, 0), "mm"))
# Save this plot
ggsave(file.path(figures.folder, "mlsa-edyears-bw.png"),
    plot = bandwidth.plot,
    units = "cm", width = fig.width, height = fig.height)


################################################################################
## Calculate the CCT bandwidth option plot, school leaving age

# Loop across the bandwidth values.
by.value <- 0.25
max.value <- 10
bw_grid <- seq(2 * by.value, max.value - 2 * by.value, by.value)
ests <- se <- numeric(length(bw_grid))
for (i in seq_along(bw_grid)) {
    print(i)
    fit <- rdrobust(
        x = full.data$birthyearmonth,
        fuzzy = full.data$agefinishededuc,
        y = log(full.data$soc_median_hourly),
        c = cutoff,
        h = bw_grid[i])
    ests[i] <- fit$coef[3]
    se[i]   <- fit$se[3]
}

# Get the CCT bandiwdth.
fit <- rdrobust(
    x = full.data$birthyearmonth,
    fuzzy = full.data$agefinishededuc,
    y = log(full.data$soc_median_hourly),
    c = cutoff)
cct.bw <- fit$bws[1,1]

# Save the output.
df_plot <- data.frame(bw = bw_grid,
    est = ests,
    lo = ests - 1.96*se,
    hi = ests + 1.96*se)
# Plot the outcome.
bandwidth.plot <- df_plot %>%
    ggplot(aes(x = bw, y = est)) +
    geom_vline(xintercept = cct.bw, linetype = "dotted") +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.25) +
    geom_point(colour = colour.list[4]) +
    geom_line(colour = colour.list[4]) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.2, fill = colour.list[4]) +
    annotate("text", x = cct.bw  * 1.05, y = 0.36,
        label = "CCT \nbandwidth",
        size = 4.25,  hjust = 0, vjust = 1) +
    theme_bw() +
    scale_x_continuous(expand = c(0, 0),
        name = "RD bandwidth choice, years to Sept 1957 birth date",
        limits = c(0, max.value - 2 * by.value),
        breaks = seq(0, max.value, by = 4 * by.value)) +
    scale_y_continuous(expand = c(0, 0),
        name = "",
        limits = c(-0.1, 0.41),
        breaks = seq(-1, 3, by = 0.1)) +
    ggtitle(TeX(r"(Fuzzy RD, Age left school on Log hourly wage)")) +
    theme(plot.title = element_text(size = rel(1), hjust = 0),
        plot.title.position = "plot",
        plot.margin = unit(c(0.5, 3, 1, 0), "mm"))
# Save this plot
ggsave(file.path(figures.folder, "mlsa-edage-bw.png"),
    plot = bandwidth.plot,
    units = "cm", width = fig.width, height = fig.height)
