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
# The standard, linear, IV estimator package.
library(ivreg)
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
    option = "none", half.line.slope = 1, half.line.intercept = 0){
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
            binscatter.ggplot <- binscatter.ggplot +
                geom_smooth(data = binscatter.data$data.dots,
                    aes(x = x, y = (half.line.intercept + fit * half.line.slope)),
                        se = FALSE, colour = "orange", size = 1, linetype = "solid")
        }
    # Return the ggplot of this.
    return(binscatter.ggplot)
}

# Load my coded version of ORIV.
source("oriv.R")


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
## Effect on education years Ed PGI -> Ed Years

# Define the set of demographic controls.
control.list <- paste0(
    "+ sex_male + factor(visityear) + poly(recruitedage, 3) + factor(sibling_count)")

causal_edpgi.reg <- function(outcome, outcome.name, input.data, control.vars){
    ## Run ech of the 6 models.
    # (1) raw OLS, Ed PGI -> outcome.
    raw.ols <- lm(formula(paste0(outcome, "~ 1 + edpgi_all_imputed_self")),
        data = input.data)
    # (2) OLS + controls, Ed PGI -> outcome
    controls.ols <- lm(formula(paste0(outcome, "~ 1 + edpgi_all_imputed_self",
        control.vars)),
        data = input.data)
    # (3) OLS (correlational ORIV), Ed PGI -> outcome
    oriv.ols <- GORIV(formula(paste0(outcome, "~ 1 ", control.vars)),
        "edpgi_all_imputed_self", "edpgi_exclude_imputed_self",
        IV.list = NULL, control.list = NULL,
        IID = "eid",
        data = input.data)
    # (4) Reduced form, random component -> outcome
    reducedform.iv <- lm(formula(paste0(outcome,
        "~ 1 + edpgi_exclude_imputed_random + edpgi_all_imputed_parental",
        control.vars)),
        data = input.data)
    # (5) Regular IV, Ed PGI -> outcome
    regular.iv <- feols(formula(paste0(
        outcome, "~ 1 + edpgi_all_imputed_parental",
        control.vars, "| edpgi_all_imputed_self ~ edpgi_exclude_imputed_random")),
        data = input.data)
    # (6) ORIV, Ed PGI -> outcome
    oriv.iv <- GORIV(formula(paste0(outcome, "~ 1 ", control.vars)),
        "edpgi_all_imputed_self", "edpgi_exclude_imputed_self",
        IV.list = c(
            "edpgi_all_imputed_random", "edpgi_exclude_imputed_random"),
        control.list = c(
            "edpgi_all_imputed_parental", "edpgi_exclude_imputed_parental"),
        IID = "eid",
        data = input.data)
    # Get the point estimates, and SEs, in a row.
    table.point <- c(
        coeftable(raw.ols)["edpgi_all_imputed_self", "Estimate"],
        coeftable(controls.ols)["edpgi_all_imputed_self", "Estimate"],
        coeftable(oriv.ols)["fit_PGI_MAIN", "Estimate"],
        coeftable(reducedform.iv)["edpgi_exclude_imputed_random", "Estimate"],
        coeftable(regular.iv)["fit_edpgi_all_imputed_self", "Estimate"],
        coeftable(oriv.iv)["fit_PGI_MAIN", "Estimate"]) %>% 
        #signif(digits.no) %>%
        style_sigfig(digits = digits.no, decimals = digits.no) %>%
        format(scientific = FALSE) %>%
        c(outcome.name, ., as.character(round(regular.iv$nobs)))
    table.se <- c(
        coeftable(raw.ols)["edpgi_all_imputed_self", "Std. Error"],
        coeftable(controls.ols)["edpgi_all_imputed_self", "Std. Error"],
        coeftable(oriv.ols)["fit_PGI_MAIN", "Std. Error"],
        coeftable(reducedform.iv)["edpgi_exclude_imputed_random", "Std. Error"],
        coeftable(regular.iv)["fit_edpgi_all_imputed_self", "Std. Error"],
        coeftable(oriv.iv)["fit_PGI_MAIN", "Std. Error"]) %>% 
        #signif(digits.no) %>%
        style_sigfig(digits = digits.no, decimals = digits.no) %>%
        format(scientific = FALSE) %>%
        paste0("(", ., ")") %>%
        c(" ", ., " ")
    # Return both rows.
    return(rbind(table.point, table.se))
}

# Calculate the relevant rows.
reg.table <- rbind(
    causal_edpgi.reg("edyears", "Education years", analysis.data, control.list),
    causal_edpgi.reg("edqual_highered", "University degree", analysis.data, control.list),
    causal_edpgi.reg("log(soc_median_hourly)", "Occupation hourly wage", analysis.data, control.list),
    causal_edpgi.reg("log(soc_median_annual)", "Occupation annual income", analysis.data, control.list)
)
# Save the LaTeX table
reg.table %>%
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

## Same, but now with other PGI as controls.
pgi_control.list <- paste0(
    "+ adhd_pgi + asthma_pgi + bipolar_pgi + bmi_pgi + height_pgi",
    "+ schizophrenia_pgi + t2diabetes_pgi",
    "+ sex_male + factor(visityear) + poly(recruitedage, 3) + factor(sibling_count)")

# Calculate the relevant rows.
pgi_controls.table <- rbind(
    causal_edpgi.reg("edyears", "Education years", analysis.data, pgi_control.list),
    causal_edpgi.reg("edqual_highered", "University degree", analysis.data, pgi_control.list),
    causal_edpgi.reg("log(soc_median_hourly)", "Occupation hourly wage", analysis.data, pgi_control.list),
    causal_edpgi.reg("log(soc_median_annual)", "Occupation annual income", analysis.data, pgi_control.list)
)

# Save the LaTeX table
pgi_controls.table %>%
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
        file = file.path(tables.folder, "genetic-effects-pgicontrols.tex"))


################################################################################
## Figure: Ed PGI -> Ed years

# Extract point-estimates from the OLS + ORIV estimates.
edyears.est <- causal_edpgi.reg("edyears", "Education years", analysis.data, control.list)
edyears.ols <- c(edyears.est[1, 2], edyears.est[2, 2]) %>%
    str_replace("\\(", "") %>%
    str_replace("\\)", "") %>%
    as.numeric()
edyears.iv <- c(edyears.est[1, 7], edyears.est[2, 7]) %>%
    str_replace("\\(", "") %>%
    str_replace("\\)", "") %>%
    as.numeric()

# Show correlation between Ed PGI and edyears in a Bin-scatter plot.
edpgi_edyears.plot <- analysis.data %>%
    binscatter.plot(data = .,
        "edpgi_all_imputed_self", "edyears", colour.list[2],
        option = "half-line",
        half.line.slope = edyears.iv[1] / edyears.ols[1],
        half.line.intercept = 2.5) +
    # Annotate OLS
    annotate("text", colour = colour.list[2],
        x = -2.5, y = 17.25,
        fontface = "bold",
        label = paste0("Raw OLS = +", edyears.ols[1], " (", edyears.ols[2], ")"),
        size = 4.25, hjust = 0, vjust = 0) +
    # Annotate IV
    annotate("text", colour = "orange",
        x = 0.5, y = 11.25,
        fontface = "bold",
        label = paste0("ORIV = +", edyears.iv[1], " (", edyears.iv[2], ")"),
        size = 4.25, hjust = 0.5, vjust = 0) +
    annotate("curve", colour = "orange",
        x = 1.75, y = 15.25,
        xend = 2, yend = 14.5,
        linewidth = 1,
        curvature = -0.25,
        arrow = arrow(length = unit(0.25, 'cm'))) +
    theme_bw() +
    scale_x_continuous(expand = c(0, 0),
        name = "Ed PGI, s.d. units",
        breaks = seq(-5, 5, by = 1),
        limits = c(-3, 3)) +
    scale_y_continuous(expand = c(0, 0.1),
        name = "",
        limits = c(11, 17.75),
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
## Figure: Ed PGI -> Ed years

# Extract point-estimates from the OLS + ORIV estimates.
earnings.est <- causal_edpgi.reg(
    "log(soc_median_hourly)", "Occupation hourly wage",
    analysis.data, control.list)
earnings.ols <- c(earnings.est[1, 2], earnings.est[2, 2]) %>%
    str_replace("\\(", "") %>%
    str_replace("\\)", "") %>%
    as.numeric()
earnings.iv <- c(earnings.est[1, 7], earnings.est[2, 7]) %>%
    str_replace("\\(", "") %>%
    str_replace("\\)", "") %>%
    as.numeric()

# Show correlation between Ed PGI and edyears in a Bin-scatter plot.
edpgi_earnings.plot <- analysis.data %>%
    mutate(log_soc_median_hourly = log(soc_median_hourly)) %>%
    binscatter.plot(data = ., "edpgi_all_imputed_self", "soc_median_hourly", colour.list[3],
        option = "half-line",
        half.line.slope = earnings.iv[1] / earnings.ols[1],
        half.line.intercept = 2) +
    # Annotate OLS
    annotate("text", colour = colour.list[3],
        x = -2.5, y = 24.25,
        fontface = "bold",
        label = paste0("Raw OLS = +", earnings.ols[1], " (", earnings.ols[2], ")"),
        size = 4.25, hjust = 0, vjust = 0) +
    # Annotate IV
    annotate("text", colour = "orange",
        x = 0.5, y = 13.75,
        fontface = "bold",
        label = paste0("ORIV = +", earnings.iv[1], " (", earnings.iv[2], ")"),
        size = 4.25, hjust = 0.5, vjust = 0) +
    annotate("curve", colour = "orange",
        x = 1.7, y = 21,
        xend = 2, yend = 19.5,
        linewidth = 1,
        curvature = -0.25,
        arrow = arrow(length = unit(0.25, 'cm'))) +
    theme_bw() +
    scale_x_continuous(expand = c(0, 0),
        name = "Ed PGI, s.d. units",
        breaks = seq(-5, 5, by = 1),
        limits = c(-3, 3)) +
    scale_y_continuous(expand = c(0, 0.1),
        name = "",
        limits = c(12.5, 25.1),
        breaks = seq(0, 50, by = 2.5)) +
    ggtitle("Occupation Hourly Wages, £") +
    theme(plot.title = element_text(size = rel(1), hjust = 0),
        plot.title.position = "plot",
        plot.margin = unit(c(0.5, 3, 0, 0), "mm"))
# Save this plot
ggsave(file.path(figures.folder, "edpgi-earnings-causal.png"),
    plot = edpgi_earnings.plot,
    units = "cm", width = fig.width, height = fig.height)


################################################################################
## Same table as before, with FEs.

outcome <- "edyears"
outcome.name <- "Education years"

# Function, adding FEs.
causal_edpgi_FEs.reg <- function(outcome, outcome.name, input.data){
    # Define the set of demographic controls.
    control.vars <- paste0(
        "+ adhd_pgi + asthma_pgi + bipolar_pgi + bmi_pgi + height_pgi",
        "+ schizophrenia_pgi + t2diabetes_pgi + sex_male",
        "+ factor(visityear) + poly(recruitedage, 3) + factor(sibling_count) + urban")
    ## Run ech of the 6 models.
    # (1) raw OLS, Ed PGI -> outcome.
    raw.ols <- feols(
        formula(paste0(outcome, "~ 1 + edpgi_all_imputed_self | famid")),
        data = input.data)
    # (2) OLS + controls, Ed PGI -> outcome
    controls.ols <- feols(
        formula(paste0(outcome, "~ 1 + edpgi_all_imputed_self",
        control.vars, " | famid")),
        data = input.data)
    # (3) OLS (correlational ORIV), Ed PGI -> outcome
    oriv.ols <- GORIV(formula(paste0(outcome, "~ 1 ", control.vars)),
        "edpgi_all_imputed_self", "edpgi_exclude_imputed_self",
        IV.list = NULL, control.list = NULL,
        IID = "eid", FID = "famid",
        data = input.data)
    # (4) Reduced form, random component -> outcome
    reducedform.iv <- feols(formula(paste0(outcome,
        "~ 1 + edpgi_exclude_imputed_random + edpgi_all_imputed_parental",
        control.vars, " | famid")),
        data = input.data)
    # (5) Regular IV, Ed PGI -> outcome
    regular.iv <- feols(formula(paste0(
        outcome, "~ 1 + edpgi_all_imputed_parental",
        control.vars, " | famid", "| edpgi_all_imputed_self ~ edpgi_exclude_imputed_random")),
        data = input.data)
    # (6) ORIV, Ed PGI -> outcome
    oriv.iv <- GORIV(formula(paste0(outcome, "~ 1 ", control.vars)),
        "edpgi_all_imputed_self", "edpgi_exclude_imputed_self",
        IV.list = c(
            "edpgi_all_imputed_random", "edpgi_exclude_imputed_random"),
        control.list = c(
            "edpgi_all_imputed_parental", "edpgi_exclude_imputed_parental"),
        IID = "eid", FID = "famid",
        data = input.data)
    # Get the point estimates, and SEs, in a row.
    table.point <- c(
        coeftable(raw.ols)["edpgi_all_imputed_self", "Estimate"],
        coeftable(controls.ols)["edpgi_all_imputed_self", "Estimate"],
        coeftable(oriv.ols)["fit_PGI_MAIN", "Estimate"],
        coeftable(reducedform.iv)["edpgi_exclude_imputed_random", "Estimate"],
        coeftable(regular.iv)["fit_edpgi_all_imputed_self", "Estimate"],
        coeftable(oriv.iv)["fit_PGI_MAIN", "Estimate"]) %>% 
        #signif(digits.no) %>%
        style_sigfig(digits = digits.no, decimals = digits.no) %>%
        format(scientific = FALSE) %>%
        c(outcome.name, ., as.character(round(regular.iv$nobs)))
    table.se <- c(
        coeftable(raw.ols)["edpgi_all_imputed_self", "Std. Error"],
        coeftable(controls.ols)["edpgi_all_imputed_self", "Std. Error"],
        coeftable(oriv.ols)["fit_PGI_MAIN", "Std. Error"],
        coeftable(reducedform.iv)["edpgi_exclude_imputed_random", "Std. Error"],
        coeftable(regular.iv)["fit_edpgi_all_imputed_self", "Std. Error"],
        coeftable(oriv.iv)["fit_PGI_MAIN", "Std. Error"]) %>% 
        #signif(digits.no) %>%
        style_sigfig(digits = digits.no, decimals = digits.no) %>%
        format(scientific = FALSE) %>%
        paste0("(", ., ")") %>%
        c(" ", ., " ")
    # Return both rows.
    return(rbind(table.point, table.se))
}

# Calculate the relevant rows.
reg_fe.table <- rbind(
    causal_edpgi_FEs.reg("edyears", "Education years", analysis.data),
    causal_edpgi_FEs.reg("edqual_highered", "University degree", analysis.data),
    causal_edpgi_FEs.reg("log(soc_median_hourly)", "Occupation hourly wage", analysis.data),
    causal_edpgi_FEs.reg("log(soc_median_annual)", "Occupation annual income", analysis.data) 
)

# Save the LaTeX table
reg_fe.table %>%
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
        file = file.path(tables.folder, "genetic-effects-FEs.tex"))


################################################################################
## Results for sample with at least one parent (i.e., more accurate imputation).

# Calculate the relevant rows.
names(analysis.data)
parental.data <- analysis.data %>%
    filter(father_present + mother_present > 0)
parental.table <- rbind(
    causal_edpgi.reg("edyears", "Education years", parental.data),
    causal_edpgi.reg("edqual_highered", "University degree", parental.data),
    causal_edpgi.reg("log(soc_median_hourly)", "Occupation hourly wage", parental.data),
    causal_edpgi.reg("log(soc_median_annual)", "Occupation annual income", parental.data) 
)

# Save the LaTeX table
parental.table %>%
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
        file = file.path(tables.folder, "genetic-effects-parental.tex"))
