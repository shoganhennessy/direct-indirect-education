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


################################################################################
## Import relevant UKB data. 

# Load the pre-cleaned UKB panel data.
ukb.data <- data.folder %>%
    file.path("ukb-cleaned-pheno.csv") %>%
    read_csv()

# Get the sibling imputed analysis sample.
analysis.data <- ukb.data %>%
    filter(analysis_sample == 1)


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
    # Generate summary data, for analysis sample (MEAN)
    summary.mean <- given.data %>%
        summarise_all(mean, na.rm = TRUE) %>%
        pivot_longer(cols = everything(),
            names_to = "variable", values_to = "mean")
    summary.sd <- given.data %>%
        summarise_all(sd, na.rm = TRUE) %>%
        pivot_longer(cols = everything(),
            names_to = "variable", values_to = "sd")
    summary.mean$sd <- summary.sd$sd
    return(summary.mean)
}

# Summarise analysis sample
analysis_summary.data <- analysis.data %>%
    transmute(
        `Male` = sex_male,
        `Age`  = recruitedage,
        `Race $=$ White` = genetic_euroancest,
        # Genetic variables.
        `Ed PGI` = edpgi_self,
        `Ed PGI, parental mean` = edpgi_parents,
        # Education variables
        `Education years` = edyears,
        # Income variables
        `Occ. hourly wage` = soc_mean_hourly,
        `Household income < 18k` = householdincome_less18k,
        `Household income 18--31k` = householdincome_18to31k,
        `Household income 31--52k` = householdincome_31to52k,
        `Household income 52--100k` = householdincome_52to100k,
        `Household income 100k <` = householdincome_above100k,
        # Designators
        `Any siblings` = sibling_present,
        `Count siblings` = sibling_count) %>%
    summary.table() %>%
    transmute(variable = variable,
        mean_analysis = mean,
        sd_analysis = sd)
# Summarise entire sample
all_summary.data <- ukb.data %>%
    # Make Ed PGI for parents missing in the entire data file summary table.
    mutate(edpgi_parents = NA) %>%
    transmute(
        `Male` = sex_male,
        `Age`  = recruitedage,
        `Race $=$ White` = genetic_euroancest,
        # Genetic variables.
        `Ed PGI` = edpgi_self,
        `Ed PGI, parental mean` = edpgi_parents,
        # Education variables
        `Education years` = edyears,
        # Income variables
        `Occ. hourly wage` = soc_mean_hourly,
        `Household income < 18k` = householdincome_less18k,
        `Household income 18--31k` = householdincome_18to31k,
        `Household income 31--52k` = householdincome_31to52k,
        `Household income 52--100k` = householdincome_52to100k,
        `Household income 100k <` = householdincome_above100k,
        # Designators
        `Any siblings` = sibling_present,
        `Count siblings` = sibling_count) %>%
    summary.table()

# Combine into one file.
summary.data <- analysis_summary.data
summary.data$mean_all <- all_summary.data$mean
summary.data$sd_all <- all_summary.data$sd

# Save the summary table as LaTeX.
summary.data %>%
    xtable(rownames = FALSE,
        out = file.path(tables.folder, "ukb-summary.tex"))

# Save the observation count, to host in files.
analysis.data %>%
    nrow() %>%
    writeLines(file.path(tables.folder, "ukb-analysis-count.txt"))
ukb.data %>%
    nrow() %>%
    writeLines(file.path(tables.folder, "ukb-total-count.txt"))


################################################################################
## Plot summaries of these data.

# Show the histogram of gene scores.
genescore.plot <- ukb.data %>%
    ggplot(aes(x = genescore_educ_euro)) +
    geom_density(aes(y = after_stat(count)), fill = colour.list[1]) +
    theme_bw() +
    geom_vline(xintercept = 0, linetype = "dashed") +
    annotate("text", colour = colour.list[1], x = -2.9, y = 2000,
        label = expression("Mean\nEd PGI = 0"),
        size = 4.25,  hjust = 0, vjust = 0) +
    scale_x_continuous(expand = c(0, 0),
        name = "Ed PGI, s.d. units",
        breaks = seq(-5, 5, by = 1)) +
    scale_y_continuous(expand = c(0, 0),
        name = "",
        #breaks = seq(0, 1, by = 0.05),
        limits = c(0, 2500)) +
    ggtitle(TeX(r"(Observations, \textit{N})")) +
    theme(plot.title = element_text(size = rel(1), hjust = 0),
        plot.title.position = "plot",
        plot.margin = unit(c(0.5, 3, 0, 0), "mm"))
# Save this plot
ggsave(file.path(figures.folder, "genescore-hist.png"),
    plot = genescore.plot,
    units = "cm", width = fig.width, height = fig.height)
# Save this plot
ggsave(file.path(figures.folder, "genescore-hist-presentation.png"),
    plot = genescore.plot,
    units = "cm", width = presentation.width, height = presentation.height)

# Show the histogram of ed years.
mean.edyears <- ukb.data %>% pull(indiv_edyears) %>% mean(na.rm = TRUE)
edyears.plot <- ukb.data %>%
    filter(indiv_edyears >= 10) %>%
    ggplot(aes(x = indiv_edyears)) +
    geom_bar(aes(y = after_stat(count)),
        fill = colour.list[2], colour = 1) +
    #geom_vline(xintercept = mean.edyears, linetype = "dashed") +
    annotate("text", colour = colour.list[2], x = 17, y = 1850,
        label = "Mean Ed years = 14 years \n(College)",
        size = 4.25, hjust = 1, vjust = 0) +
    theme_bw() +
    scale_x_continuous(expand = c(0, 0.01),
        name = "Education Years",
        breaks = seq(0, 20, by = 1)) +
    scale_y_continuous(expand = c(0, 0),
        name = "",
        breaks = seq(0, 3000, by = 250),
        limits = c(0, 2240)) +
    ggtitle(TeX(r"(Observations, \textit{N})")) +
    theme(plot.title = element_text(size = rel(1), hjust = 0),
        plot.title.position = "plot",
        plot.margin = unit(c(0.5, 3, 0, 0), "mm"))
# Save this plot
ggsave(file.path(figures.folder, "edyears-hist.png"),
    plot = edyears.plot,
    units = "cm", width = fig.width, height = fig.height)

# Show the histogram of annual income.
mean.earnings <- ukb.data %>% pull(indiv_earnings_real) %>% mean(na.rm = TRUE)
earnings.plot <- ukb.data %>%
    ggplot(aes(x = indiv_earnings_real / 1000)) +
    geom_density(aes(y = after_stat(count)), fill = colour.list[3]) +
    geom_vline(xintercept = mean.earnings / 1000, linetype = "dashed") +
    annotate("text", colour = colour.list[3], x = 100, y = 50,
        label = "Mean Individual \nEarnings = $80,000",
        size = 4.25,  hjust = 0, vjust = 0) +
    theme_bw() +
    scale_x_continuous(expand = c(0, 0.5),
        name = "Annual Earnings, $ thousands",
        limits = c(0, 300),
        breaks = seq(0, 300, by = 50)) +
    scale_y_continuous(expand = c(0, 0),
        name = "",
        breaks = seq(0, 100, by = 10),
        limits = c(0, 62.5)) +
    ggtitle(TeX(r"(Observations, \textit{N})")) +
    theme(plot.title = element_text(size = rel(1), hjust = 0),
        plot.title.position = "plot",
        plot.margin = unit(c(0.5, 3, 0, 0), "mm"))
# Save this plot
ggsave(file.path(figures.folder, "earnings-hist.png"),
    plot = earnings.plot,
    units = "cm", width = fig.width, height = fig.height)

# Show correlation between edyears and income
ukb.data %>%
    filter(10 <= indiv_edyears, indiv_edyears <= 17) %>%
    lm(log(indiv_earnings_real) ~ 1 + indiv_edyears, data = .) %>%
    summary() %>%
    print()
# Scatter plot of edyears + earnings.
edyears_earnings.plot <- ukb.data %>%
    filter(indiv_edyears >= 10) %>%
    group_by(indiv_edyears) %>%
    summarise(indiv_earnings_real = mean(indiv_earnings_real, na.rm = TRUE)) %>%
    ggplot(aes(x = indiv_edyears, y = indiv_earnings_real / 1000)) +
    geom_point() +
    geom_smooth(method = "loess",
        colour = colour.list[2], fill = colour.list[2], size = 2) +
    annotate("text", colour = colour.list[2],
        x = 13.5, y = 137.5,
        fontface = "bold",
        label = "Slope = +11.5%",
        size = 4.25, hjust = 0.5, vjust = 0) +
    theme_bw() +
    scale_x_continuous(
        expand = c(0, 0),
        name = "Education Years",
        breaks = seq(0, 20, by = 1)) +
    scale_y_continuous(expand = c(0, 1),
        name = "",
        limits = c(0, 150),
        breaks = seq(0, 300, by = 25)) +
    ggtitle("Annual Earnings, $ thousands") +
    theme(plot.title = element_text(size = rel(1), hjust = 0),
        plot.title.position = "plot",
        plot.margin = unit(c(0.5, 3, 0, 0), "mm"))
# Save this plot
ggsave(file.path(figures.folder, "edyears-earnings.png"),
    plot = edyears_earnings.plot,
    units = "cm", width = fig.width, height = fig.height)
# Save this plot
ggsave(file.path(figures.folder, "presentation-edyears-earnings.png"),
    plot = edyears_earnings.plot,
    units = "cm", width = presentation.width, height = 0.75 * fig.height)

## Show correlation between edyears and income
# Get mean eductaion gain between high-school and uni degree.
ukb.data %>%
    filter(indiv_edyears >= 12, degree_reported >= 2) %>%
    mutate(college_degree = as.integer(degree_reported >= 5)) %>%
    lm(log(indiv_earnings_real) ~ 1 + college_degree, data = .) %>%
    summary() %>%
    print()
# SHow this in a plot.
college_earnings.plot <- ukb.data %>%
    filter(indiv_edyears >= 12, degree_reported >= 2) %>%
    mutate(college_degree = as.character(degree_reported >= 5)) %>%
    group_by(college_degree) %>%
    summarise(indiv_earnings_real = mean(indiv_earnings_real, na.rm = TRUE)) %>%
    ungroup() %>%
    ggplot(aes(x = college_degree, y = indiv_earnings_real / 1000)) +
    geom_bar(stat = "identity", colour = 1, fill = colour.list[3]) +
    annotate("text", colour = colour.list[3],
        x = 0.5, y = 100,
        label = expression("Mean ≈ +45%"),
        size = 4.25, hjust = 0, vjust = 0) +
    theme_bw() +
    scale_x_discrete(
        #expand = c(0.05, 0.05),
        name = "Higher Education?",
        labels = c("TRUE" = "Undergrad+", "FALSE" = "Sec. School+")) +
    scale_y_continuous(expand = c(0, 0),
        name = "",
        limits = c(0, 150),
        breaks = seq(0, 150, by = 25)) +
    ggtitle("Annual Earnings, $ thousands") +
    theme(plot.title = element_text(size = rel(1), hjust = 0),
        plot.title.position = "plot",
        plot.margin = unit(c(0.5, 3, 0, 0), "mm"))
# Save this plot
ggsave(file.path(figures.folder, "college-earnings.png"),
    plot = college_earnings.plot,
    units = "cm", width = fig.width, height = fig.height)

# Show correlation between gene score and edyears
# First-stage OLS no controls 0.66519
ukb.data %>%
    filter(10 <= indiv_edyears, indiv_edyears <= 17) %>%
    lm(indiv_edyears ~ 1 + genescore_educ_euro, data = .) %>%
    summary() %>%
    print()
# Show in a Bin-scatter plot.
genescore_edyears.plot <- ukb.data %>%
    filter(10 <= indiv_edyears, indiv_edyears <= 17) %>%
    binscatter.plot(data = ., "genescore_educ_euro", "indiv_edyears",
        colour.list[2]) +
    annotate("text", colour = colour.list[2],
        x = 0.875, y = 10.5,
        fontface = "bold",
        label = ("Slope = +0.61 (0.02) Ed years"),
        size = 4.25, hjust = 0.5, vjust = 0) +
    theme_bw() +
    scale_x_continuous(expand = c(0, 0),
        name = "Ed PGI, s.d. units",
        breaks = seq(-5, 5, by = 1),
        limits = c(-3, 3)) +
    scale_y_continuous(expand = c(0, 0.1),
        name = "",
        limits = c(10, 17),
        breaks = seq(0, 20, by = 1)) +
    ggtitle("Education Years") +
    theme(plot.title = element_text(size = rel(1), hjust = 0),
        plot.title.position = "plot",
        plot.margin = unit(c(0.5, 3, 0, 0), "mm"))
# Save this plot
ggsave(file.path(figures.folder, "genescore-edyears.png"),
    plot = genescore_edyears.plot,
    units = "cm", width = fig.width, height = fig.height)
# Save this plot
ggsave(file.path(figures.folder, "genescore-edyears-wider.png"),
    plot = genescore_edyears.plot,
    units = "cm", width = presentation.width, height = 0.85 * presentation.height)

# Generate a version with the (implied) causal design.
genescore_edyears_causal.plot <- ukb.data %>%
    filter(10 <= indiv_edyears, indiv_edyears <= 17) %>%
    binscatter.plot(data = ., "genescore_educ_euro", "indiv_edyears",
        colour.list[2], option = "half-line") +
    annotate("text", colour = colour.list[2],
        x = 0.875, y = 10.5,
        fontface = "bold",
        label = ("Slope = +0.61 (0.02) Ed years"),
        size = 4.25, hjust = 0.5, vjust = 0) +
    # Add on the part with a lower slope.
    annotate("text", colour = "orange",
        x = 1.125, y = 11.25,
        fontface = "bold",
        label = ("Slope = +0.30 (Young et al., 2022)"),
        size = 4.25, hjust = 0.5, vjust = 0) +
    annotate("curve", colour = "orange",
        x = 1.75, y = 14.75,
        xend = 2, yend = 14.00,
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
        limits = c(10, 17),
        breaks = seq(0, 20, by = 1)) +
    ggtitle("Education Years") +
    theme(plot.title = element_text(size = rel(1), hjust = 0),
        plot.title.position = "plot",
        plot.margin = unit(c(0.5, 3, 0, 0), "mm"))
# Save this plot.
ggsave(file.path(figures.folder, "genescore-edyears-causal.png"),
    plot = genescore_edyears_causal.plot,
    units = "cm", width = presentation.width, height = 0.85 * presentation.height)

# Show correlation between gene score and income
# Reduced-form OLS no controls 0.100560
# First-stage OLS no controls 0.66519
ukb.data %>%
    filter(indiv_earnings_real <= 250000) %>%
    lm(log(indiv_earnings_real) ~ 1 + genescore_educ_euro, data = .) %>%
    summary() %>%
    print()
# Reduced-form OLS with controls 0.056507
genescore_earnings.plot <- ukb.data %>%
    mutate(indiv_earnings_real = indiv_earnings_real / 1000) %>%
    binscatter.plot(data = ., "genescore_educ_euro", "indiv_earnings_real",
        colour.list[3]) +
    annotate("text", colour = colour.list[3],
        x = 1.5, y = 15,
        fontface = "bold",
        label = ("Slope = +8.8% (0.008)"),
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
    ggtitle("Annual Earnings, $ thousands") +
    theme(plot.title = element_text(size = rel(1), hjust = 0),
        plot.title.position = "plot",
        plot.margin = unit(c(0.5, 3, 0, 0), "mm"))
# Save this plot
ggsave(file.path(figures.folder, "genescore-earnings.png"),
    plot = genescore_earnings.plot,
    units = "cm", width = fig.width, height = fig.height)
# Save this plot
ggsave(file.path(figures.folder, "genescore-earnings-wider.png"),
    plot = genescore_earnings.plot,
    units = "cm", width = presentation.width, height = 0.85 * presentation.height)

# Generate a version with the (implied) causal design.
genescore_earnings_causal.plot <- ukb.data %>%
    mutate(indiv_earnings_real = indiv_earnings_real / 1000) %>%
    binscatter.plot(data = ., "genescore_educ_euro", "indiv_earnings_real",
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
    ggtitle("Annual Earnings, $ thousands") +
    theme(plot.title = element_text(size = rel(1), hjust = 0),
        plot.title.position = "plot",
        plot.margin = unit(c(0.5, 3, 0, 0), "mm"))
# Save this plot
ggsave(file.path(figures.folder, "genescore-earnings-causal.png"),
    plot = genescore_earnings_causal.plot,
    units = "cm", width = presentation.width, height = 0.85 * presentation.height)


# Correlation between gene score and cognitive measures later in life
# Example of a direct effect.
lm(log(cogtot) ~ 1 + genescore_educ_euro, data = ukb.data) %>% summary()
genescore_cogtot.plot <- ukb.data %>%
    binscatter.plot(data = ., "genescore_educ_euro", "cogtot", "orange") +
    annotate("text", colour = "orange",
        x = 0, y = 37.5,
        label = expression("Slope ≈ +2.9%"),
        size = 4.25,  hjust = 0, vjust = 0) +
    theme_bw() +
    scale_x_continuous(expand = c(0, 0),
        name = "Ed PGI, s.d. units",
        breaks = seq(-5, 5, by = 1)) +
    scale_y_continuous(expand = c(0, 1),
        name = "",
        limits = c(0, 40),
        breaks = seq(0, 40, by = 10)) +
    ggtitle("Performance on a Cognitive Exam, at age 55--65") +
    theme(plot.title = element_text(size = rel(1), hjust = 0),
        plot.title.position = "plot",
        plot.margin = unit(c(0.5, 3, 0, 0), "mm"))
# Save this plot
ggsave(file.path(figures.folder, "genescore-cogtot.png"),
    plot = genescore_cogtot.plot,
    units = "cm", width = fig.width, height = fig.height)
# Save this plot
ggsave(file.path(figures.folder, "genescore-cogtot-wider.png"),
    plot = genescore_cogtot.plot,
    units = "cm", width = presentation.width, height = fig.height)

# Compare education years.
ukb.data %>%
    dplyr::select(indiv_edyears, degree_reported) %>%
    table(exclude = NULL)

# See dist of edyears.
ukb.data %>%
    pull(indiv_edyears) %>%
    table(exclude = NULL)


################################################################################
## Correlation matrix, including Ed PGI

# Show how Ed PGI is related to demographics, lm(Z ~ 1 + X)
demographic.reg <- ukb.data %>%
    transmute(
        genescore_educ_euro = genescore_educ_euro,
        `Childhood SES: Family poor`            = family_poor,
        `Childhood SES: Family moved`           = family_move,
        `Childhood SES: Family financial help`  = family_finhelp,
        `Childhood SES: Father missing`         = father_missing,
        `Childhood SES: Father unemployed`      = father_unemp,
        `Childhood SES: Father manual labourer` = father_manualjob,
        `Childhood SES: Parents smoked`         = parents_smoke,
        `Childhood SES: Childhood head injury`  = child_headinjury,
        `Age`                                   = indiv_agey,
        `No. children`                          = child,
        `Size household`                        = hhres,
        `Father Ed years`                       = ifelse(is.na(father_edyears), 0, father_edyears),
        `Mother Ed years`                       = ifelse(is.na(mother_edyears), 0, mother_edyears),
        `Mean Parent Ed years`                  = parent_edyears,
        `Female`                                = gender_female) %>%
    lm(reformulate(".", response = "genescore_educ_euro"), data = .)
print(summary(demographic.reg))
# Plot the coefficients.
demographic.plot <- modelsummary::modelplot(demographic.reg,
        coef_omit = "Intercept", colour = "blue", size = 1) +
    theme_bw() +
    geom_vline(xintercept = 0, linetype = "dashed") +
    scale_x_continuous(expand = c(0.005, 0.005),
        name = TeX(r"(Association Estimates, ED PGI$_i = \beta\,' \, X_i + \epsilon_i$)"),
        breaks = seq(-1, 1, by = 0.05)) +
    ggtitle(TeX("Demographic Information")) +
    theme(plot.title = element_text(size = rel(1), hjust = 0),
        plot.title.position = "plot",
        plot.margin = unit(c(0.5, 3, 0, 0), "mm"),
        axis.text.y = element_text(hjust = 0))
demographic.plot
# Save this plot
ggsave(file.path(figures.folder, "demographic-plot.png"),
    plot = demographic.plot,
    units = "cm",
    width = 1.5 * presentation.width, height = 1.25 * presentation.height)

# Calculate the correlation plot.
genescore_correlation.matrix <- ukb.data %>%
    filter(!is.na(genescore_educ_euro)) %>%
    # Select the gene scores
    transmute(
        EA           = genescore_educ_euro,
        Cognition    = genescore_gencog_euro,
        BMI          = genescore_bmi_euro,
        Alcoholism   = genescore_alcohol_euro,
        Alzheimers   = genescore_alzh_euro,
        Wellbeing    = genescore_wellbeing_euro,
        Bipolar      = genescore_bipolar_euro,
        ADHD         = genescore_adhd_euro,
        Longevity    = genescore_longevity_euro,
        Antisocial   = genescore_antisocial_euro,
        Depression   = genescore_mdepressived_euro,
        Anxiety      = genescore_anxietyfactor_euro) %>%
    cor()
# Plot the matrix, only the lower left.
genescore_correlation.matrix[
    upper.tri(genescore_correlation.matrix, diag = FALSE)] <- NA
genescore_correlation.plot <- genescore_correlation.matrix %>%
    ggcorrplot::ggcorrplot(hc.order = FALSE,
        type = "full")
# Save this file
ggsave(file.path(figures.folder, "genescore-correlation.png"),
    plot = genescore_correlation.plot,
    units = "cm", width = 1.5 * fig.width, height = 1.25 * fig.height)


################################################################################
## Plot comparison of causal med (1) OLS (2) OLS + controls (3) DML.

# Causal medation package, Imai Keele Yamamoto (2010)
library(mediation)

# Function to extract the coefficients from the mediate objects
effects.extract <- function(mediate.reg, model.name){
    # Compile the mediation regression results.
    reg.summary <- summary(mediate.reg)
    # Get the total effect estimates.
    total.est <- reg.summary$tau.coef
    total.ci.upper <- as.numeric(reg.summary$tau.ci[1])
    total.ci.lower <- as.numeric(reg.summary$tau.ci[2])
    # Get the direct effect estimates.
    direct.est <- reg.summary$z.avg
    direct.ci.upper <- as.numeric(reg.summary$z.avg.ci[1])
    direct.ci.lower <- as.numeric(reg.summary$z.avg.ci[2])
    # Get the indirect effect estimates.
    indirect.est <- reg.summary$d.avg
    indirect.ci.upper <- as.numeric(reg.summary$d.avg.ci[1])
    indirect.ci.lower <- as.numeric(reg.summary$d.avg.ci[2])
    # Get the percent mediated estimates.
    permediated.est <- reg.summary$n.avg
    permediated.ci.upper <- as.numeric(reg.summary$n.avg.ci[1])
    permediated.ci.lower <- as.numeric(reg.summary$n.avg.ci[2])
    # Put it all into a dataframe
    data.return <- data.frame(
        effect = c("Total", "Direct", "Indirect", "Percent Mediated"),
        pointest = c(total.est, direct.est, indirect.est, permediated.est),
        upperest = c(total.ci.upper, direct.ci.upper, indirect.ci.upper, permediated.ci.upper),
        lowerest = c(total.ci.lower, direct.ci.lower, indirect.ci.lower, permediated.ci.lower))
    # Label it with the model name
    data.return <- data.return %>% mutate(model = model.name)
    return(data.return)
}

# Put relevant data into matrix form
Z <- ukb.data %>% pull(genescore_educ_euro)
D <- ukb.data %>% pull(indiv_edyears)
Y <- ukb.data %>% transmute(Y = log(indiv_earnings_real)) %>% pull(Y)
X <- ukb.data %>%
    dplyr::select(c(
        # Regular demographics
        survey_year, indiv_agey, parent_edyears, gender_female,
        # Childhood SES
        family_poor, family_move, family_finhelp,
        father_missing, father_unemp, father_manualjob,
        parents_smoke, child_headinjury)) %>%
    as.matrix()
# Consider D = higher education (binary needed)
#Z <- as.integer(Z >= 0)
D <- as.integer(D >= 14)

# Define number of samples to bootstrap over.
boot.samples <- 10^4

# Define the simple OLS model
mediation_ols.reg <- 
    mediate(
        lm(D ~ (1 + Z)),
        lm(Y ~ (1 + Z * D)),
        treat = "Z", mediator = "D",
        robustSE = FALSE, sims = boot.samples)
print(summary(mediation_ols.reg))

# Define the OLS model with linear controls
mediation_controls.reg <- 
    mediate(
        lm(D ~ (1 + Z) + X),
        lm(Y ~ (1 + Z * D) + X),
        treat = "Z", mediator = "D",
        robustSE = FALSE, sims = boot.samples)
print(summary(mediation_controls.reg))

# Define the OLS model with non-linear controls
mediation_nonparametric.reg <-
    mediate(
        lm(D ~ (1 + Z) * X),
        lm(Y ~ (1 + Z * D) * X),
        treat = "Z", mediator = "D",
        robustSE = FALSE, sims = boot.samples)
print(summary(mediation_nonparametric.reg))

# Get all the results together.
effects.data <- rbind(
    effects.extract(mediation_ols.reg, "OLS"),
    effects.extract(mediation_controls.reg, "+ Controls"),
    effects.extract(mediation_nonparametric.reg, "+ Non-Linear"))

# Plot in a bar chart.
mediation.plot <- effects.data %>%
    filter(effect != "Percent Mediated") %>%
    ggplot(aes(
        fill = factor(effect, levels = c("Total", "Direct", "Indirect")),
        x = model)) +
    geom_bar(aes(y = pointest),
        stat = "identity", position = "dodge", colour = "black") +
    geom_errorbar(aes(ymin = lowerest, ymax = upperest),
        size = 2 / 3,
        stat = "identity", position = position_dodge(0.9), width = 1 / 3) +
    theme_bw() +
    geom_hline(yintercept = 0, alpha = 0.5) +
    scale_fill_discrete(
        name = "",
        limits = c("Total", "Direct", "Indirect"),
        type = colour.list[c(3, 1, 2)]) +
    scale_x_discrete(
        name = "",
        limits = c("OLS", "+ Controls", "+ Non-Linear")) +
    scale_y_continuous(expand = c(0, 0, 0.01, 0),
        #limits = c(0, 0.12),
        name = "",
        breaks = seq(0, 0.5, by = 0.025)) +
    ggtitle("Estimate, Percent on Earnings") +
    theme(plot.title = element_text(size = rel(1), hjust = 0),
        plot.title.position = "plot",
        plot.margin = unit(c(0.5, 3, 0.25, 0), "mm"),
        legend.position = c(0.775, 0.925),
        legend.direction = "horizontal")
# Save this file
ggsave(file.path(figures.folder, "mediation-plot.png"),
    plot = mediation.plot,
    units = "cm",
    width = presentation.width * 1.125, height = presentation.height * 1.125)

# Plot in a bar chart, just the total effect
mediation_total.plot <- effects.data %>%
    #filter(effect == "Total") %>%
    mutate(
        pointest = ifelse(effect == "Total", pointest, 0),
        upperest = ifelse(effect == "Total", upperest, 0),
        lowerest = ifelse(effect == "Total", lowerest, 0)) %>%
    ggplot(aes(
        fill = factor(effect, levels = c("Total", "Direct", "Indirect")),
        x = model)) +
    geom_bar(aes(y = pointest),
        stat = "identity", position = "dodge", colour = "black") +
    geom_errorbar(aes(ymin = lowerest, ymax = upperest),
        size = 2 / 3,
        stat = "identity", position = position_dodge(0.9), width = 1 / 3) +
    theme_bw() +
    geom_hline(yintercept = 0, alpha = 0.5) +
    scale_fill_discrete(
        name = "",
        limits = c("Total", "Direct", "Indirect"),
        type = colour.list[c(3, 1, 2)]) +
    scale_x_discrete(
        name = "",
        limits = c("OLS", "+ Controls", "+ Non-Linear")) +
    scale_y_continuous(expand = c(0, 0, 0.01, 0),
        #limits = c(0, 0.12),
        name = "",
        breaks = seq(0, 0.5, by = 0.025)) +
    ggtitle("Estimate, Percent on Earnings") +
    theme(plot.title = element_text(size = rel(1), hjust = 0),
        plot.title.position = "plot",
        plot.margin = unit(c(0.5, 3, 0.25, 0), "mm"),
        legend.position = c(0.775, 0.925),
        legend.direction = "horizontal")
# Save this file
ggsave(file.path(figures.folder, "mediation-total-plot.png"),
    plot = mediation_total.plot,
    units = "cm",
    width = presentation.width * 1.125, height = presentation.height * 1.125)

# Plot in a bar chart, just the percent mediatied.
mediation_percent.plot <- effects.data %>%
    filter(effect == "Percent Mediated") %>%
    ggplot(aes(x = model)) +
    geom_bar(aes(y = pointest), fill = "orange",
        stat = "identity", position = "dodge", colour = "black") +
    geom_errorbar(aes(ymin = lowerest, ymax = upperest),
        size = 2 / 3,
        stat = "identity", position = position_dodge(0.9), width = 1 / 3) +
    theme_bw() +
    geom_hline(yintercept = 0, alpha = 0.5) +
    scale_x_discrete(
        name = "",
        limits = c("OLS", "+ Controls", "+ Non-Linear")) +
    scale_y_continuous(expand = c(0, 0, 0.01, 0),
        limits = c(0, 1),
        name = "",
        breaks = seq(0, 1, by = 0.1)) +
    ggtitle("Estimate, Percent Mediated Through Higher Education") +
    theme(plot.title = element_text(size = rel(1), hjust = 0),
        plot.title.position = "plot",
        plot.margin = unit(c(0.5, 3, 0.25, 0), "mm"),
        legend.position = c(0.775, 0.925),
        legend.direction = "horizontal")
# Save this file
ggsave(file.path(figures.folder, "mediation-percent-plot.png"),
    plot = mediation_percent.plot,
    units = "cm",
    width = presentation.width * 1.125, height = presentation.height * 1.125)


################################################################################
## Selection model estimates of causal mediation

# Point estimate with the Heckman parametric selection model.
library(sampleSelection)
library(boot)
boot.samples <- 10^4
# Put relevant data into matrix form
Z <- ukb.data %>% pull(genescore_educ_euro)
D <- ukb.data %>% pull(indiv_edyears)
Y <- ukb.data %>% transmute(Y = log(indiv_earnings_real)) %>% pull(Y)
X <- ukb.data %>%
    dplyr::select(c(
        survey_year, indiv_agey, parent_edyears, gender_female)) %>%
    as.matrix()
# Consider D = higher education (binary needed)
#Z <- as.integer(Z >= 0)
D <- as.integer(D >= 14)

# Define a function that outputs a parametric selection model's estimates of
# (1) total effect (2) First-stage (3) direct (4) indirect
selection_mediation.reg <- function(data, type = "parametric"){
    # Get the relevant columns from the provided data.
    Z <- data[["Z"]]
    D <- data[["D"]]
    Y <- data[["Y"]]
    X <- as.matrix(data[4:dim(data)[2]])
    # Unified first-stage
    firststage.reg <- probit(D ~ (1 + Z) * X)
    # Unified total effect
    totaleffect.reg <- lm(Y ~ (1 + Z) * X)
    ## Get the prediction models by type
    if (type == "parametric"){
        # Get the second-stage by using a parametric control function.
        # Unobserved part estimated with a Heckman (switching) selection approach.
        selection.reg <- selection(
            # First-stage
            as.formula(D ~ (1 + Z) * X),
            # Second-stage
            list(as.formula(Y ~ (1 + Z)), as.formula(Y ~ (1 + Z) * X)),
                method = "2step")
        # Generate predictions
        EY_Zz_D0 <- predict(
            selection.reg, newdata = data.frame(Z = Z, D = 0, X = X))[, 1]
        EY_Zz_D1 <- predict(
            selection.reg, newdata = data.frame(Z = Z, D = 1, X = X))[, 2]
        EY_Z1_Dd <- predict(
            selection.reg, newdata = data.frame(Z = Z + 0.5, D = D, X = X))
        EY_Z1_Dd <- ifelse(D == 1, EY_Z1_Dd[, 2], EY_Z1_Dd[, 1])
        EY_Z0_Dd <- predict(
            selection.reg, newdata = data.frame(Z = Z - 0.5, D = D, X = X))
        EY_Z0_Dd <- ifelse(D == 1, EY_Z0_Dd[, 2], EY_Z0_Dd[, 1])
    }
    else if (type == "semi-parametric"){
        controlfun <- (D - predict(firststage.reg, type = "response"))
        selection.reg <- lm(Y ~ (1 + Z * D) + controlfun * D)
        # Generate predictions
        EY_Zz_D0 <- predict(
            selection.reg, newdata = data.frame(Z = Z, D = 0, X = X,
                controlfun = controlfun), rankdeficient = "simple")
        EY_Zz_D1 <- predict(
            selection.reg, newdata = data.frame(Z = Z, D = 1, X = X,
                controlfun = controlfun), rankdeficient = "simple")
        EY_Z1_Dd <- predict(
            selection.reg, newdata = data.frame(Z = Z + 0.5, D = D, X = X,
                controlfun = controlfun), rankdeficient = "simple")
        EY_Z0_Dd <- predict(
            selection.reg, newdata = data.frame(Z = Z - 0.5, D = D, X = X,
                controlfun = controlfun), rankdeficient = "simple")
    }
    # Get the total effect by prediction
    totaleffect.est <- mean(
        (predict(totaleffect.reg, newdata = data.frame(Z = Z + 0.5, X = X)) -
            predict(totaleffect.reg, newdata = data.frame(Z = Z - 0.5, X = X))))
    # Get the first-stage by prediction
    complier.score <- (predict(firststage.reg, 
        newdata = data.frame(Z = Z + 0.5, X = X), type = "response") -
            predict(firststage.reg,
                newdata = data.frame(Z = Z - 0.5, X = X), type = "response"))
    firststage.est <- mean(complier.score)
    # Direct Effect, E[ Y_i(1, D) - Y_i(0, D)].
    direct.gains <- EY_Z1_Dd - EY_Z0_Dd
    direct.est <- mean(direct.gains)
    # Indirect effect, E[ Y_i(Z, 1) - Y_i(Z, 0) \times (D_i(1) - D_i(0)) ]
    treatment.gains <- EY_Zz_D1 - EY_Zz_D0
    indirect.est <- mean(complier.score * treatment.gains)
    treatment.est <- weighted.mean(treatment.gains, complier.score)
    # Percent mediated
    # note: probably does not follow a normal at boundary -> bootstrap inapprop.
    mediated.est <- indirect.est / totaleffect.est
    # Return the relevant statistics
    return(c(totaleffect.est, firststage.est,
        direct.est, indirect.est, mediated.est, treatment.est))
}

# Show how the function works
selection.example <- selection_mediation.reg(
    data.frame(Z = Z, D = D, Y = Y, X = X), type = "parametric")
print(selection.example)
# Or the semi-para est
selection.example <- selection_mediation.reg(
    data.frame(Z = Z, D = D, Y = Y, X = X), type = "semi-parametric")
print(selection.example)

# Define the bootstrap function.
selection_mediation.boot <- function(data, indicies){
    estimates.list <- selection_mediation.reg(
        data[indicies, ], type = "parametric")
    return(estimates.list)
}

# Define the bootstrap function.
selection_semiparametric.boot <- function(data, indicies){
    estimates.list <- selection_mediation.reg(
        data[indicies, ], type = "semi-parametric")
    return(estimates.list)
}

# Bootstrap the selection mediated estimates.
selection_mediation.bootest <-
    boot(data = data.frame(Z = Z, D = D, Y = Y, X = X),
        statistic = selection_mediation.boot, R = boot.samples)
selection_semiparametric.bootest <-
    boot(data = data.frame(Z = Z, D = D, Y = Y, X = X),
        statistic = selection_semiparametric.boot, R = boot.samples)

# Function to extract the coefficients from the bootstrapped estimates
bootstrap.extract <- function(model.bootest, model.name){
    # Extract the columns of the boot object by name
    total.boot       <- model.bootest$t[, 1]
    firststage.boot  <- model.bootest$t[, 2]
    direct.boot      <- model.bootest$t[, 3]
    indirect.boot    <- model.bootest$t[, 4]
    permediated.boot <- model.bootest$t[, 5]
    # Get the total effect estimates.
    total.est      <- mean(total.boot)
    total.se       <- sd(total.boot)
    total.ci.upper <- as.numeric(quantile(total.boot, probs = 0.975))
    total.ci.lower <- as.numeric(quantile(total.boot, probs = 0.025))
    # Get the direct effect estimates.
    direct.est      <- mean(direct.boot)
    direct.se       <- sd(direct.boot)
    direct.ci.upper <- as.numeric(quantile(direct.boot, probs = 0.975))
    direct.ci.lower <- as.numeric(quantile(direct.boot, probs = 0.025))
    # Get the indirect effect estimates.
    indirect.est      <- mean(indirect.boot)
    indirect.se       <- sd(indirect.boot)
    indirect.ci.upper <- as.numeric(quantile(indirect.boot, probs = 0.975))
    indirect.ci.lower <- as.numeric(quantile(indirect.boot, probs = 0.025))
    # Get the percent mediated estimates.
    permediated.est      <- mean(permediated.boot)
    permediated.se       <- sd(permediated.boot)
    permediated.ci.upper <- as.numeric(quantile(permediated.boot, probs = 0.975))
    permediated.ci.lower <- as.numeric(quantile(permediated.boot, probs = 0.025))
    # Put it all into a dataframe
    data.return <- data.frame(
        effect = c("Total", "Direct", "Indirect", "Percent Mediated"),
        pointest = c(total.est, direct.est, indirect.est, permediated.est),
        upperest = c(total.ci.upper, direct.ci.upper, indirect.ci.upper, permediated.ci.upper),
        lowerest = c(total.ci.lower, direct.ci.lower, indirect.ci.lower, permediated.ci.lower))
    # Label it with the model name
    data.return <- data.return %>% mutate(model = model.name)
    return(data.return)
}

# Get the regular mediation estimates of the binary model.
mediation_binary.reg <-
    mediate(
        lm(D ~ (1 + Z) * X),
        lm(Y ~ (1 + Z * D) * X),
        treat = "Z", mediator = "D",
        robustSE = FALSE, sims = boot.samples)
print(summary(mediation_binary.reg))

# Append the estimates to those extracted earlier
selection_est.data <- rbind(
    effects.extract(mediation_binary.reg, "OLS + Controls"),
    bootstrap.extract(selection_mediation.bootest,
        "Parametric Heckman\nSelection Model"),
    bootstrap.extract(selection_semiparametric.bootest,
        "Semi-Parametric\nSelection Model"))

# Plot in a bar chart.
mediation_selection.plot <- selection_est.data %>%
    filter(effect != "Percent Mediated") %>%
    ggplot(aes(
        fill = factor(effect, levels = c("Total", "Direct", "Indirect")),
        x = model)) +
    geom_bar(aes(y = pointest),
        stat = "identity", position = "dodge", colour = "black") +
    geom_errorbar(aes(ymin = lowerest, ymax = upperest),
        size = 2 / 3,
        stat = "identity", position = position_dodge(0.9), width = 1 / 3) +
    theme_bw() +
    geom_hline(yintercept = 0, alpha = 0.5) +
    scale_fill_discrete(
        name = "",
        limits = c("Total", "Direct", "Indirect"),
        type = colour.list[c(3, 1, 2)]) +
    scale_x_discrete(
        name = "",
        limits = c("OLS + Controls",
            "Parametric Heckman\nSelection Model",
            "Semi-Parametric\nSelection Model")) +
    scale_y_continuous(expand = c(0, 0, 0.01, 0),
        #limits = c(-0.03, 0.175),
        name = "",
        breaks = seq(-0.5, 0.5, by = 0.025)) +
    ggtitle("Estimate, Percent on Earnings") +
    theme(plot.title = element_text(size = rel(1), hjust = 0),
        plot.title.position = "plot",
        plot.margin = unit(c(0, 3, 0.25, 0), "mm"),
        legend.position = c(0.775, 0.925),
        legend.direction = "horizontal")
# Save this file
ggsave(file.path(figures.folder, "mediation-selection-plot.png"),
    plot = mediation_selection.plot,
    units = "cm",
    width = presentation.width * 1.125, height = presentation.height * 1.125)

# Plot in a bar chart, just the percent mediatied.
mediation_selection_percent.plot <- selection_est.data %>%
    #mutate(
    #    pointest = ifelse(pointest > 1, 1, pointest),
    #    upperest = ifelse(upperest > 1, 1, upperest),
    #    lowerest = ifelse(lowerest > 1, 1, lowerest)) %>%
    filter(effect == "Percent Mediated") %>%
    ggplot(aes(x = model)) +
    geom_bar(aes(y = pointest), fill = "orange",
        stat = "identity", position = "dodge", colour = "black") +
    geom_errorbar(aes(ymin = lowerest, ymax = upperest),
        size = 2 / 3,
        stat = "identity", position = position_dodge(0.9), width = 1 / 3) +
    theme_bw() +
    geom_hline(yintercept = 0, alpha = 0.5) +
    geom_hline(yintercept = 1, alpha = 0.5, linetype = "dotted") +
    scale_x_discrete(
        name = "",
        limits = c("OLS + Controls",
            "Parametric Heckman\nSelection Model",
            "Semi-Parametric\nSelection Model")) +
    scale_y_continuous(expand = c(0, 0, 0.05, 0),
        #limits = c(0, 1),
        name = "",
        breaks = seq(-2, 2, by = 0.1)) +
    ggtitle("Estimate, Percent Mediated Through Higher Education") +
    theme(plot.title = element_text(size = rel(1), hjust = 0),
        plot.title.position = "plot",
        plot.margin = unit(c(0.5, 3, 0.25, 0), "mm"),
        legend.position = c(0.775, 0.925),
        legend.direction = "horizontal")
# Save this file
ggsave(file.path(figures.folder, "mediation-selection-percent-plot.png"),
    plot = mediation_selection_percent.plot,
    units = "cm",
    width = presentation.width * 1.125, height = presentation.height * 1.125)

# Question: What are the implied returns to education estimates?
