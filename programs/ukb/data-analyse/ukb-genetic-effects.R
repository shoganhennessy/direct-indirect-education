#!/usr/bin/R
## Senan Hogan-Hennessy, 24 March 2025
## UKB data -> Genetic effects, Ed PGI random component -> Ed Years + Income.
set.seed(47)
print(Sys.time())

# Functions for data manipulation and visualisation
library(tidyverse)
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
binscatter.plot <- function(data, x, y, colour.name, option = "none"){
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
            binscatter.ggplot <- binscatter.ggplot +
                geom_smooth(data = binscatter.data$data.dots,
                aes(x = x, y = (6 + fit / 2)), se = FALSE,
                colour = "orange", size = 1, linetype = "solid")
        }
    # Return the ggplot of this.
    return(binscatter.ggplot)
}

# Load my coded version of ORIV.
source("oriv.R")

# Define  function to take the name of a variable, and estimate every relevant model.


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

outcome <- "edyears"
outcome.name <- "Education years"

causal_edpgi.reg <- function(outcome, outcome.name){
    # Define the set of demographic controls.
    control.vars <- paste0(
        "+ adhd_pgi + asthma_pgi + bipolar_pgi + bmi_pgi + height_pgi",
        "+ schizophrenia_pgi + t2diabetes_pgi + sex_male",
        "+ factor(visityear) + poly(recruitedage, 3) + factor(sibling_count) + urban")
    ## Run ech of the 6 models.
    # (1) raw OLS, Ed PGI -> outcome.
    raw.ols <- lm(formula(paste0(outcome, "~ 1 + edpgi_all_imputed_self")),
        data = analysis.data)
    # (2) OLS + controls, Ed PGI -> outcome
    controls.ols <- lm(formula(paste0(outcome, "~ 1 + edpgi_all_imputed_self",
        control.vars)),
        data = analysis.data)
    # (3) OLS (correlational ORIV), Ed PGI -> outcome
    oriv.ols <- GORIV(formula(paste0(outcome, "~ 1 ", control.vars)),
        "edpgi_all_imputed_self", "edpgi_exclude_imputed_self",
        IV.list = NULL, control.list = NULL,
        IID = "eid",
        data = analysis.data)
    # (4) Reduced form, random component -> outcome
    reducedform.iv <- lm(formula(paste0(outcome,
        "~ 1 + edpgi_all_imputed_random + edpgi_all_imputed_parental",
        control.vars)),
        data = analysis.data)
    # (5) Regular IV, Ed PGI -> outcome
    regular.iv <- feols(formula(paste0(
        outcome, "~ 1 + edpgi_all_imputed_parental",
        control.vars, "| edpgi_all_imputed_self ~ edpgi_exclude_imputed_random")),
        data = analysis.data)
    summary(regular.iv)
    # (6) ORIV, Ed PGI -> outcome
    oriv.iv <- GORIV(formula(paste0(outcome, "~ 1 ", control.vars)),
        "edpgi_all_imputed_self", "edpgi_exclude_imputed_self",
        IV.list = c(
            "edpgi_all_imputed_random", "edpgi_exclude_imputed_random"),
        control.list = c(
            "edpgi_all_imputed_parental", "edpgi_exclude_imputed_parental"),
        IID = "eid",
        data = analysis.data)
    # Get the point estimates, and SEs, in a row.
    table.point <- c(
        coeftable(raw.ols)["edpgi_all_imputed_self", "Estimate"],
        coeftable(controls.ols)["edpgi_all_imputed_self", "Estimate"],
        coeftable(oriv.ols)["fit_PGI_MAIN", "Estimate"],
        coeftable(reducedform.iv)["edpgi_all_imputed_random", "Estimate"],
        coeftable(regular.iv)["fit_edpgi_all_imputed_self", "Estimate"],
        coeftable(oriv.iv)["fit_PGI_MAIN", "Estimate"]) %>% 
        signif(digits.no) %>%
        format(scientific = FALSE) %>%
        c(outcome.name, ., as.character(round(regular.iv$nobs)))
    table.se <- c(
        coeftable(raw.ols)["edpgi_all_imputed_self", "Std. Error"],
        coeftable(controls.ols)["edpgi_all_imputed_self", "Std. Error"],
        coeftable(oriv.ols)["fit_PGI_MAIN", "Std. Error"],
        coeftable(reducedform.iv)["edpgi_all_imputed_random", "Std. Error"],
        coeftable(regular.iv)["fit_edpgi_all_imputed_self", "Std. Error"],
        coeftable(oriv.iv)["fit_PGI_MAIN", "Std. Error"]) %>% 
        signif(digits.no) %>%
        format(scientific = FALSE) %>%
        paste0("(", ., ")") %>%
        c(" ", ., " ")
    # Return both rows.
    return(rbind(table.point, table.se))
}


# Calculate the relevant rows.
reg.table <- rbind(
    causal_edpgi.reg("edyears", "Education years"),
    causal_edpgi.reg("edqual_highered", "University degree"),
    causal_edpgi.reg("log(soc_mean_hourly)", "Occupation hourly wage"),
    causal_edpgi.reg("log(soc_mean_annual)", "Occupation annual income")
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


################################################################################
## Figure: Ed PGI -> Ed years

#TODO: extract point-estimates from the OLS + ORIV estimates.
# Show correlation between Ed PGI and edyears
analysis.data %>%
    lm(edyears ~ 1 + edpgi_all_imputed_self, data = .) %>%
    summary() %>%
    print()
# Show in a Bin-scatter plot.
edpgi_edyears.plot <- analysis.data %>%
    binscatter.plot(data = ., "edpgi_all_imputed_self", "edyears", colour.list[1]) +
    annotate("text", colour = colour.list[1],
        x = -2.5, y = 16,
        fontface = "bold",
        label = ("Slope = +0.93 (0.02) \n Ed years"),
        size = 4.25, hjust = 0, vjust = 0) +
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
ggsave(file.path(figures.folder, "edpgi-edyears.png"),
    plot = edpgi_edyears.plot,
    units = "cm", width = fig.width, height = fig.height)