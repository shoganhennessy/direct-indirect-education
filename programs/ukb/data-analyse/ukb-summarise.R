#!/usr/bin/R
## Senan Hogan-Hennessy, 24 March 2025
## Summarise UKB data.
set.seed(47)
print(Sys.time())

# Functions for data manipulation and visualisation
library(tidyverse)
# Functions for tables into TeX
library(xtable)
# Library for better colour choice.
library(ggthemes)
# Library for equations in plots
library(latex2exp)
# Library for Bin Scatter
library(binsreg)
# Define number of digits in tables and graphs
digits.no <- 3
# Size for figures
fig.width <- 8.5
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
presentation.folder <- file.path("..", "..", "..", "presentation",
    "presentation-files", "figures")
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


################################################################################
## Import relevant UKB data. 

# Load the pre-cleaned UKB panel data.
ukb.data <- data.folder %>%
    file.path("ukb-cleaned-pheno.csv") %>%
    read_csv() 

# Get the sibling imputed analysis sample.
analysis.data <- ukb.data %>%
    filter(analysis_sample == 1)

#! Test.
analysis.data %>% head(100) %>% print()
analysis.data %>% names() %>% print()

# Show me the first row of the dataset, to describe in the data section.
analysis.data %>%
    head(1) %>%
    select(visityear, sex_male, recruitedage, jobcode_soc,
        soc_median_hourly, soc_median_annual, hours_workweek) %>%
    print()


################################################################################
## Show characteristics of the samples, in a summary table.

# WANT: Summary table of this flavour:
#        | UKB sibling sample | UKB entire sample
#  ----------------------------------------------
#  Var 1 | Mean               | Mean             
#        | (SD)               | (SD)             
#  ...                                            
#  Var N | Mean               | Mean             
#        | (SD)               | (SD)             
#        |                    |                  
#        | Obs                | Obs.

summary.table <- function(given.data){
    # Subset for useful variables.
    table.data <- given.data %>%
        transmute(
            # Basic demographics.
            "\\\\[-1.8ex] \\textit{Demographics:}" = NA,
            "Male" = sex_male,
            "Age"  = recruitedage,
            "Lives in city" = urban,
            "Ethnicity $=$ European" = genetic_euroancest,
            "Any siblings also in UKB?"     = sibling_present,
            "Count of siblings in UKB" = sibling_count,
            # Genetic variables.
            "\\\\[-1.8ex]\\hline \\\\[-1.8ex] \\textit{Genetic Measures:}" = NA,
            "Ed PGI"                        = edpgi_all_imputed_self,
            "Ed PGI, imputed parental mean" = edpgi_all_imputed_parental,
            "Other PGI: ADHD"               = adhd_pgi,
            "Other PGI: Asthma"             = asthma_pgi,
            "Other PGI: Bipolar"            = bipolar_pgi,
            "Other PGI: BMI"                = bmi_pgi,
            "Other PGI: Diabetes type 2"    = t2diabetes_pgi,
            "Other PGI: Height"             = height_pgi,
            "Other PGI: Schizophrenia"      = schizophrenia_pgi,
            # Education variables
            "\\\\[-1.8ex]\\hline \\\\[-1.8ex] \\textit{Education:}" = NA,
            "Education years" = edyears,
            "Age left education" = agefinishededuc,
            "Qualification, University degree" = edqual_highered,
            "Qualification, A-Levels" = edqual_alevels,
            "Qualification, GCSEs" = edqual_gcses,
            "Qualification, Professional degree" = edqual_professional,
            "Qualification, Vocational degree" = edqual_vocational,
            "Qualification, No official qualifications" = edqual_minimum,
            # Income variables
            "\\\\[-1.8ex]\\hline \\\\[-1.8ex] \\textit{Labour Market Outcomes:}" = NA,
            "Occupation hourly wage, \\pounds"             = soc_median_hourly,
            "Occupation annual income, thousands \\pounds" = soc_median_annual,
            "Average hours worked, per week"               = hours_workweek,
            "Household income, $< \\pounds 18k$"           = householdincome_less18k,
            "Household income, $\\pounds 18-31k$"          = householdincome_18to31k,
            "Household income, $\\pounds 31-52k$"          = householdincome_31to52k,
            "Household income, $\\pounds 52-100k$"         = householdincome_52to100k,
            "Household income, $\\pounds 100k <$"          = householdincome_above100k)
    # Generate summary data, for provided data (MEAN)
    summary.mean <- table.data %>%
        summarise_all(mean, na.rm = TRUE) %>%
        pivot_longer(cols = everything(),
            names_to = "variable", values_to = "mean")
    summary.sd <- table.data %>%
        summarise_all(sd, na.rm = TRUE) %>%
        pivot_longer(cols = everything(),
            names_to = "variable", values_to = "sd")
    summary.mean$sd <- summary.sd$sd
    return(summary.mean)
}

# Summarise analysis sample
analysis_summary.data <- analysis.data %>%
    summary.table() %>%
    transmute(variable = variable,
        mean_analysis = mean,
        sd_analysis = sd)
# Summarise entire sample
all_summary.data <- ukb.data %>%
    # Use Ed PGI for entire sample in table, and missing for parents
    mutate(
        edpgi_all_imputed_self = edpgi_all_imputed_self,
        edpgi_all_imputed_parental = NA) %>%
    summary.table()

# Combine into one file.
summary.data <- analysis_summary.data
summary.data$mean_all <- all_summary.data$mean
summary.data$sd_all <- all_summary.data$sd

# Save the summary table as LaTeX.
summary.data %>%
    xtable() %>%
    print(
        sanitize.colnames.function = identity,
        sanitize.text.function = identity,
        NA.string = " ",
        include.colnames = FALSE,
        include.rownames = FALSE,
        only.contents = TRUE,
        hline.after = NULL,
        format.args = list(big.mark = ","),
        file = file.path(tables.folder, "ukb-summary.tex"))

# Save the observation count, to host in files.
analysis.data %>%
    nrow() %>%
    prettyNum(big.mark = ",", scientific = FALSE) %>%
    writeLines(file.path(tables.folder, "ukb-analysis-count.txt"))
ukb.data %>%
    nrow() %>%
    prettyNum(big.mark = ",", scientific = FALSE) %>%
    writeLines(file.path(tables.folder, "ukb-total-count.txt"))


################################################################################
## Plot data summaries of these data.

# Show the histogram of Ed PGIs.
edpgi.plot <- analysis.data %>%
    ggplot(aes(x = edpgi_all_imputed_self)) +
    geom_density(aes(y = after_stat(count)), fill = colour.list[1]) +
    theme_bw() +
    geom_vline(xintercept = 0, linetype = "dashed") +
    annotate("text", colour = colour.list[1],
        x = -3.55, y = 8250,
        fontface = "bold",
        label = "Mean Ed PGI \n               = 0",
        size = 4.25,  hjust = 0, vjust = 0) +
    scale_x_continuous(expand = c(0, 0),
        name = "Ed PGI, s.d. units",
        breaks = seq(-5, 5, by = 1)) +
    scale_y_continuous(expand = c(0, 0),
        name = "",
        limits = c(0, 10010)) +
    ggtitle(TeX(r"(Observations)")) +
    theme(plot.title = element_text(size = rel(1), hjust = 0),
        plot.title.position = "plot",
        plot.margin = unit(c(0.5, 3, 0, 0), "mm"))
# Save this plot
ggsave(file.path(figures.folder, "edpgi-hist.png"),
    plot = edpgi.plot,
    units = "cm", width = fig.width, height = fig.height)
# Save this plot
ggsave(file.path(presentation.folder, "edpgi-hist.png"),
    plot = edpgi.plot,
    units = "cm", width = presentation.width, height = presentation.height)

# Show the histogram of ed years.
mean.edyears <- analysis.data %>% pull(edyears) %>% mean(na.rm = TRUE)
analysis.data %>% pull(edyears) %>% table(exclude = NULL) %>% print()
edyears.plot <- analysis.data %>%
    ggplot(aes(x = edyears)) +
    geom_bar(aes(y = after_stat(count)),
        fill = colour.list[2], colour = 1) +
    #geom_vline(xintercept = mean.edyears, linetype = "dashed") +
    annotate("text", colour = colour.list[2], x = 16, y = 7500,
        fontface = "bold",
        label = paste0("Mean Ed years \n= ", round(mean.edyears, 2),
            " (Sec. school)"),
        size = 4.25, hjust = 1, vjust = 0) +
    theme_bw() +
    scale_x_continuous(expand = c(0.01, 0.01),
        limits = c(8.5, 18.5),
        name = "Education Years",
        breaks = seq(0, 20, by = 1)) +
    scale_y_continuous(expand = c(0, 0),
        name = "",
        breaks = seq(0, 10000, by = 1000),
        limits = c(0, 9000)) +
    ggtitle(TeX(r"(Observations, \textit{N})")) +
    theme(plot.title = element_text(size = rel(1), hjust = 0),
        plot.title.position = "plot",
        plot.margin = unit(c(0.5, 3, 0, 0), "mm"))
# Save this plot
ggsave(file.path(figures.folder, "edyears-hist.png"),
    plot = edyears.plot,
    units = "cm", width = fig.width, height = fig.height)

# Show the histogram of annual income.
mean.earnings <- analysis.data %>% pull(soc_median_annual) %>% mean(na.rm = TRUE)
earnings.plot <- analysis.data %>%
    ggplot(aes(x = soc_median_annual)) +
    geom_density(aes(y = after_stat(count)), fill = colour.list[3]) +
    geom_vline(xintercept = mean.earnings, linetype = "dashed") +
    #annotate("text", colour = colour.list[3], x = 100, y = 50,
    #    label = "Mean Individual \nEarnings = $80,000",
    #    size = 4.25,  hjust = 0, vjust = 0) +
    theme_bw() +
    scale_x_continuous(expand = c(0, 0.5),
        name = "Annual Income, $ thousands",
        #limits = c(0, 300),
        breaks = seq(0, 300, by = 50)) +
    scale_y_continuous(expand = c(0, 0),
        name = "") +
        #breaks = seq(0, 100, by = 10),
        #limits = c(0, 62.5)
    ggtitle(TeX(r"(Observations, \textit{N})")) +
    theme(plot.title = element_text(size = rel(1), hjust = 0),
        plot.title.position = "plot",
        plot.margin = unit(c(0.5, 3, 0, 0), "mm"))
# Save this plot
ggsave(file.path(figures.folder, "earnings-hist.png"),
    plot = earnings.plot,
    units = "cm", width = fig.width, height = fig.height)

# Show correlation between edyears and income
edyears.reg <- lm(log(soc_median_annual) ~ 1 + edyears, data = analysis.data)
edyears.point <- coef(summary(edyears.reg))["edyears", "Estimate"]
edyears.se <- coef(summary(edyears.reg))["edyears", "Std. Error"]
# Scatter plot of edyears + earnings.
edyears_earnings.plot <- analysis.data %>%
    group_by(edyears) %>%
    summarise(soc_median_annual = mean(soc_median_annual, na.rm = TRUE)) %>%
    ggplot(aes(x = edyears, y = soc_median_annual)) +
    geom_smooth(method = "loess",
        colour = colour.list[2], fill = colour.list[2], size = 2) +
    geom_point(colour = "black", size = 2) +
    annotate("text", colour = colour.list[2],
        x = 13.5, y = 40,
        fontface = "bold",
        label = paste0("Slope = +", round(100 * edyears.point, 2),
            "% \n ", "(", round(edyears.se, 3), ")"),
        size = 4.25, hjust = 1, vjust = 0) +
    theme_bw() +
    scale_x_continuous(
        expand = c(0.01, 0.01),
        limits = c(8, 18),
        name = "Education Years",
        breaks = seq(0, 20, by = 1)) +
    scale_y_continuous(expand = c(0, 1),
        name = "",
        limits = c(0, 50),
        breaks = seq(0, 100, by = 10)) +
    ggtitle("Annual Income, £ thousands") +
    theme(plot.title = element_text(size = rel(1), hjust = 0),
        plot.title.position = "plot",
        plot.margin = unit(c(0.5, 3, 0, 0), "mm"))
# Save this plot
ggsave(file.path(figures.folder, "edyears-earnings.png"),
    plot = edyears_earnings.plot,
    units = "cm", width = fig.width, height = fig.height)
# Save this plot
ggsave(file.path(presentation.folder, "edyears-earnings.png"),
    plot = edyears_earnings.plot,
    units = "cm", width = presentation.width, height = presentation.height)

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
# Save this plot
ggsave(file.path(presentation.folder, "edpgi-edyears.png"),
    plot = edpgi_edyears.plot,
    units = "cm", width = presentation.width, height = presentation.height)

# Show correlation between Ed PGI and income
# Reduced-form OLS no controls 0.100560
# First-stage OLS no controls 0.66519
analysis.data %>%
    lm(log(soc_median_annual) ~ 1 + edpgi_all_imputed_self, data = .) %>%
    summary() %>%
    print()
# Reduced-form OLS with controls 0.056507
edpgi_earnings.plot <- analysis.data %>%
    mutate(log_soc_median_annual = log(soc_median_annual)) %>%
    binscatter.plot(data = ., "edpgi_all_imputed_self", "log_soc_median_annual",
        colour.list[3]) +
    annotate("text", colour = colour.list[3],
        x = 1.5, y = 15,
        fontface = "bold",
        label = ("Slope = +8.1% (0.004)"),
        size = 4.25, hjust = 0.5, vjust = 0) +
    theme_bw() +
    scale_x_continuous(expand = c(0, 0),
        name = "Ed PGI, s.d. units",
        limits = c(-3, 3),
        breaks = seq(-5, 5, by = 1)) +
    scale_y_continuous(expand = c(0, 1),
        name = "",
        limits = c(0, 150),
        breaks = seq(0, 300, by = 25)) +
    ggtitle("Annual Income, £ thousands") +
    theme(plot.title = element_text(size = rel(1), hjust = 0),
        plot.title.position = "plot",
        plot.margin = unit(c(0.5, 3, 0, 0), "mm"))
# Save this plot
ggsave(file.path(figures.folder, "edpgi-earnings.png"),
    plot = edpgi_earnings.plot,
    units = "cm", width = fig.width, height = fig.height)
# Save this plot
ggsave(file.path(presentation.folder, "edpgi-earnings.png"),
    plot = edpgi_earnings.plot,
    units = "cm", width = presentation.width, height = presentation.height)

# Generate a version with the (implied) causal design.
edpgi_earnings_causal.plot <- analysis.data %>%
    mutate(soc_median_annual = soc_median_annual) %>%
    binscatter.plot(data = ., "edpgi_all_imputed_self", "soc_median_annual",
        colour.list[3]) +
    annotate("text", colour = colour.list[3],
        x = 1.5, y = 15,
        fontface = "bold",
        label = ("Slope = +8.8% (0.008)"),  
        size = 4.25, hjust = 0.5, vjust = 0) +
    # Add on the part with a lower slope.
    annotate("text", colour = "orange",
        x = 0.875, y = 40,
        fontface = "bold",
        label = ("Slope = ?"),
        size = 4.25, hjust = 0.5, vjust = 0) +
    annotate("curve", colour = "orange",
        x = 0.5, y = 75,
        xend = 1.1, yend = 47.5,
        linewidth = 1,
        curvature = -0.25,
        arrow = arrow(length = unit(0.25, 'cm'))) +
    theme_bw() +
    scale_x_continuous(expand = c(0, 0),
        name = "Ed PGI, s.d. units",
        limits = c(-3, 3),
        breaks = seq(-5, 5, by = 1)) +
    scale_y_continuous(expand = c(0, 1),
        name = "",
        limits = c(0, 150),
        breaks = seq(0, 300, by = 25)) +
    ggtitle("Annual Income, £ thousands") +
    theme(plot.title = element_text(size = rel(1), hjust = 0),
        plot.title.position = "plot",
        plot.margin = unit(c(0.5, 3, 0, 0), "mm"))
# Save this plot
ggsave(file.path(figures.folder, "edpgi-earnings-causal.png"),
    plot = edpgi_earnings_causal.plot,
    units = "cm", width = presentation.width, height = presentation.height)


################################################################################
## Correlation matrix, including Ed PGI

## Calculate the correlation plot.
#edpgi_correlation.matrix <- analysis.data %>%
#    # Select the Ed PGIs
#    transmute(
#        EA           = edpgi_all_imputed_self,
#        Cognition    = edpgi_gencog_euro,
#        BMI          = edpgi_bmi_euro,
#        Alcoholism   = edpgi_alcohol_euro,
#        Alzheimers   = edpgi_alzh_euro,
#        Wellbeing    = edpgi_wellbeing_euro,
#        Bipolar      = edpgi_bipolar_euro,
#        ADHD         = edpgi_adhd_euro,
#        Longevity    = edpgi_longevity_euro,
#        Antisocial   = edpgi_antisocial_euro,
#        Depression   = edpgi_mdepressived_euro,
#        Anxiety      = edpgi_anxietyfactor_euro) %>%
#    cor()
## Plot the matrix, only the lower left.
#edpgi_correlation.matrix[
#    upper.tri(edpgi_correlation.matrix, diag = FALSE)] <- NA
#edpgi_correlation.plot <- edpgi_correlation.matrix %>%
#    ggcorrplot::ggcorrplot(hc.order = FALSE,
#        type = "full")
## Save this file
#ggsave(file.path(figures.folder, "edpgi-correlation.png"),
#    plot = edpgi_correlation.plot,
#    units = "cm", width = 1.5 * fig.width, height = 1.25 * fig.height)
#