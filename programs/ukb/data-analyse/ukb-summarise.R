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

#! Test.
analysis.data %>% head(100) %>% View()


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
    # Generate summary data, for provided data (MEAN)
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
        "Male" = sex_male,
        "Age"  = recruitedage,
        "Ethnicity $=$ European" = genetic_euroancest,
        # Genetic variables.
        "Ed PGI" = edpgi_self,
        "Ed PGI, parental mean" = edpgi_parents,
        # Education variables
        "Education years" = edyears,
        # Income variables
        "Occupation hourly wage"    = soc_mean_hourly,
        "Occupation annual income"  = soc_mean_annual,
        "Average hours worked week" = hours_workweek,
        "Household income < 18k"    = householdincome_less18k,
        "Household income 18--31k"  = householdincome_18to31k,
        "Household income 31--52k"  = householdincome_31to52k,
        "Household income 52--100k" = householdincome_52to100k,
        "Household income 100k <"   = householdincome_above100k,
        # Family designators
        "Any siblings?"  = sibling_present,
        "Count of siblings" = sibling_count) %>%
    summary.table() %>%
    transmute(variable = variable,
        mean_analysis = mean,
        sd_analysis = sd)
# Summarise entire sample
all_summary.data <- ukb.data %>%
    # Make Ed PGI for parents missing in the entire data file summary table.
    mutate(edpgi_parents = NA) %>%
    transmute(
        "Male" = sex_male,
        "Age"  = recruitedage,
        "Ethnicity $=$ European" = genetic_euroancest,
        # Genetic variables.
        "Ed PGI" = edpgi_all,
        "Ed PGI, parental mean" = edpgi_parents,
        # Education variables
        "Education years" = edyears,
        # Income variables
        "Occ. hourly wage"          = soc_mean_hourly,
        "Household income $<$ 18k"  = householdincome_less18k,
        "Household income 18--31k"  = householdincome_18to31k,
        "Household income 31--52k"  = householdincome_31to52k,
        "Household income 52--100k" = householdincome_52to100k,
        "Household income 100k $<$" = householdincome_above100k,
        # Designators
        "Any siblings" = sibling_present,
        "Count siblings" = sibling_count) %>%
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

#TODO: funtionalise the above, and then make 4 panels.
analysis.data %>%
    transmute(
        # Panel A. Genetic measures.
        "Ed PGI" = edpgi_self,
        "Ed PGI, parental mean" = edpgi_parents,
        #TODO:
        # Panel B. Education.
        "Education years" = edyears,
        #TODO: "No education qualification" = as.integer(edyears >= 10)
        #TODO: "GCSEs $+$"                  = as.integer(edyears >= 10)
        #TODO: "A Levels $+$"               = as.integer(edyears >= 13)
        #TODO: "Professional degree"        = as.integer(edyears == 15)
        #TODO: "Vocational degree"          = as.integer(edyears == 19)
        #TODO: "University degree"          = as.integer(edyears == 20)
        # Panel C. Income and labour market outcomes
        "Occupation hourly wage"    = soc_mean_hourly,
        #TODO -> "Occupation annual income"    = soc_mean_hourly,
        "Household income < 18k"    = householdincome_less18k,
        "Household income 18--31k"  = householdincome_18to31k,
        "Household income 31--52k"  = householdincome_31to52k,
        "Household income 52--100k" = householdincome_52to100k,
        "Household income 100k <"   = householdincome_above100k,
        # Panel D. Demographics
        "Male"                    = sex_male,
        "Age"                     = recruitedage,
        "Ethnicity $=$ Caucasian" = genetic_euroancest,
        "Any siblings?"           = sibling_present,
        "Count of siblings"       = sibling_count)


################################################################################
## Plot data summaries of these data.

# Show the histogram of Ed PGIs.
edpgi.plot <- analysis.data %>%
    ggplot(aes(x = edpgi_self)) +
    geom_density(aes(y = after_stat(count)), fill = colour.list[1]) +
    theme_bw() +
    geom_vline(xintercept = 0, linetype = "dashed") +
    annotate("text", colour = colour.list[1],
        x = -2.9, y = 8000,
        label = TeX(r"(Mean Ed PGI $=$ 0)"),
        size = 4.25,  hjust = 0, vjust = 0) +
    scale_x_continuous(expand = c(0, 0),
        name = "Ed PGI, s.d. units",
        breaks = seq(-5, 5, by = 1)) +
    scale_y_continuous(expand = c(0, 0),
        name = "",
        limits = c(0, 10010)) +
    ggtitle(TeX(r"(Observations, \textit{N})")) +
    theme(plot.title = element_text(size = rel(1), hjust = 0),
        plot.title.position = "plot",
        plot.margin = unit(c(0.5, 3, 0, 0), "mm"))
# Save this plot
ggsave(file.path(figures.folder, "edpgi-hist.png"),
    plot = edpgi.plot,
    units = "cm", width = fig.width, height = fig.height)
# Save this plot
ggsave(file.path(figures.folder, "edpgi-hist-presentation.png"),
    plot = edpgi.plot,
    units = "cm", width = presentation.width, height = presentation.height)

# Show the histogram of ed years.
mean.edyears <- analysis.data %>% pull(edyears) %>% mean(na.rm = TRUE)
edyears.plot <- analysis.data %>%
    #filter(edyears >= 10) %>%
    ggplot(aes(x = edyears)) +
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
        name = "") +#,
        #breaks = seq(0, 3000, by = 250),
        #limits = c(0, 2240)) +
    ggtitle(TeX(r"(Observations, \textit{N})")) +
    theme(plot.title = element_text(size = rel(1), hjust = 0),
        plot.title.position = "plot",
        plot.margin = unit(c(0.5, 3, 0, 0), "mm"))
# Save this plot
ggsave(file.path(figures.folder, "edyears-hist.png"),
    plot = edyears.plot,
    units = "cm", width = fig.width, height = fig.height)

# Show the histogram of annual income.
mean.earnings <- analysis.data %>% pull(soc_mean_hourly) %>% mean(na.rm = TRUE)
earnings.plot <- analysis.data %>%
    ggplot(aes(x = soc_mean_hourly)) +
    geom_density(aes(y = after_stat(count)), fill = colour.list[3]) +
    geom_vline(xintercept = mean.earnings, linetype = "dashed") +
    #annotate("text", colour = colour.list[3], x = 100, y = 50,
    #    label = "Mean Individual \nEarnings = $80,000",
    #    size = 4.25,  hjust = 0, vjust = 0) +
    theme_bw() +
    scale_x_continuous(expand = c(0, 0.5),
        name = "Annual Earnings, $ thousands",
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
analysis.data %>%
    filter(10 <= edyears, edyears <= 17) %>%
    lm(log(soc_mean_hourly) ~ 1 + edyears, data = .) %>%
    summary() %>%
    print()
# Scatter plot of edyears + earnings.
edyears_earnings.plot <- analysis.data %>%
    filter(edyears >= 10) %>%
    group_by(edyears) %>%
    summarise(soc_mean_hourly = mean(soc_mean_hourly, na.rm = TRUE)) %>%
    ggplot(aes(x = edyears, y = soc_mean_hourly)) +
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
analysis.data %>%
    filter(edyears >= 12) %>%
    mutate(college_degree = as.integer(degree_reported >= 5)) %>%
    lm(log(soc_mean_hourly) ~ 1 + college_degree, data = .) %>%
    summary() %>%
    print()
# SHow this in a plot.
college_earnings.plot <- analysis.data %>%
    filter(edyears >= 12) %>%
    mutate(college_degree = as.character(degree_reported >= 5)) %>%
    group_by(college_degree) %>%
    summarise(soc_mean_hourly = mean(soc_mean_hourly, na.rm = TRUE)) %>%
    ungroup() %>%
    ggplot(aes(x = college_degree, y = soc_mean_hourly)) +
    geom_bar(stat = "identity", colour = 1, fill = colour.list[3]) +
    annotate("text", colour = colour.list[3],
        x = 0.5, y = 100,
        label = TeX(r"(Mean ≈ +45%)"),
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

# Show correlation between Ed PGI and edyears
# First-stage OLS no controls 0.66519
analysis.data %>%
    filter(10 <= edyears, edyears <= 17) %>%
    lm(edyears ~ 1 + edpgi_self, data = .) %>%
    summary() %>%
    print()
# Show in a Bin-scatter plot.
edpgi_edyears.plot <- analysis.data %>%
    filter(10 <= edyears, edyears <= 17) %>%
    binscatter.plot(data = ., "edpgi_self", "edyears",
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
ggsave(file.path(figures.folder, "edpgi-edyears.png"),
    plot = edpgi_edyears.plot,
    units = "cm", width = fig.width, height = fig.height)
# Save this plot
ggsave(file.path(figures.folder, "edpgi-edyears-wider.png"),
    plot = edpgi_edyears.plot,
    units = "cm", width = presentation.width, height = 0.85 * presentation.height)

# Generate a version with the (implied) causal design.
edpgi_edyears_causal.plot <- analysis.data %>%
    filter(10 <= edyears, edyears <= 17) %>%
    binscatter.plot(data = ., "edpgi_self", "edyears",
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
ggsave(file.path(figures.folder, "edpgi-edyears-causal.png"),
    plot = edpgi_edyears_causal.plot,
    units = "cm", width = presentation.width, height = 0.85 * presentation.height)

# Show correlation between Ed PGI and income
# Reduced-form OLS no controls 0.100560
# First-stage OLS no controls 0.66519
analysis.data %>%
    filter(soc_mean_hourly <= 250000) %>%
    lm(log(soc_mean_hourly) ~ 1 + edpgi_self, data = .) %>%
    summary() %>%
    print()
# Reduced-form OLS with controls 0.056507
edpgi_earnings.plot <- analysis.data %>%
    mutate(soc_mean_hourly = soc_mean_hourly) %>%
    binscatter.plot(data = ., "edpgi_self", "soc_mean_hourly",
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
ggsave(file.path(figures.folder, "edpgi-earnings.png"),
    plot = edpgi_earnings.plot,
    units = "cm", width = fig.width, height = fig.height)
# Save this plot
ggsave(file.path(figures.folder, "edpgi-earnings-wider.png"),
    plot = edpgi_earnings.plot,
    units = "cm", width = presentation.width, height = 0.85 * presentation.height)

# Generate a version with the (implied) causal design.
edpgi_earnings_causal.plot <- analysis.data %>%
    mutate(soc_mean_hourly = soc_mean_hourly) %>%
    binscatter.plot(data = ., "edpgi_self", "soc_mean_hourly",
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
ggsave(file.path(figures.folder, "edpgi-earnings-causal.png"),
    plot = edpgi_earnings_causal.plot,
    units = "cm", width = presentation.width, height = 0.85 * presentation.height)


# Correlation between Ed PGI and cognitive measures later in life
# Example of a direct effect.
lm(log(cogtot) ~ 1 + edpgi_self, data = ukb.data) %>% summary()
edpgi_cogtot.plot <- analysis.data %>%
    binscatter.plot(data = ., "edpgi_self", "cogtot", "orange") +
    annotate("text", colour = "orange",
        x = 0, y = 37.5,
        label = TeX(r"(Slope ≈ +2.9%)"),
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
ggsave(file.path(figures.folder, "edpgi-cogtot.png"),
    plot = edpgi_cogtot.plot,
    units = "cm", width = fig.width, height = fig.height)
# Save this plot
ggsave(file.path(figures.folder, "edpgi-cogtot-wider.png"),
    plot = edpgi_cogtot.plot,
    units = "cm", width = presentation.width, height = fig.height)

# Compare education years.
analysis.data %>%
    dplyr::select(edyears, degree_reported) %>%
    table(exclude = NULL)

# See dist of edyears.
analysis.data %>%
    pull(edyears) %>%
    table(exclude = NULL)


################################################################################
## Correlation matrix, including Ed PGI

# Show how Ed PGI is related to demographics, lm(Z ~ 1 + X)
demographic.reg <- analysis.data %>%
    transmute(
        edpgi_self = edpgi_self,
        "Childhood SES: Family poor"            = family_poor,
        "Childhood SES: Family moved"           = family_move,
        "Childhood SES: Family financial help"  = family_finhelp,
        "Childhood SES: Father missing"         = father_missing,
        "Childhood SES: Father unemployed"      = father_unemp,
        "Childhood SES: Father manual labourer" = father_manualjob,
        "Childhood SES: Parents smoked"         = parents_smoke,
        "Childhood SES: Childhood head injury"  = child_headinjury,
        "Age"                                   = indiv_agey,
        "No. children"                          = child,
        "Size household"                        = hhres,
        "Father Ed years"                       = ifelse(is.na(father_edyears), 0, father_edyears),
        "Mother Ed years"                       = ifelse(is.na(mother_edyears), 0, mother_edyears),
        "Mean Parent Ed years"                  = parent_edyears,
        "Female"                                = gender_female) %>%
    lm(reformulate(".", response = "edpgi_self"), data = .)
print(summary(demographic.reg))
# Plot the coefficients.
demographic.plot <- modelsummary::modelplot(demographic.reg,
        coef_omit = "Intercept", colour = "blue", size = 1) +
    theme_bw() +
    geom_vline(xintercept = 0, linetype = "dashed") +
    scale_x_continuous(expand = c(0.005, 0.005),
        name = TeX(r"(Association Estimates, Ed PGI$_i = \beta\,' \, X_i + \epsilon_i$)"),
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
edpgi_correlation.matrix <- analysis.data %>%
    # Select the Ed PGIs
    transmute(
        EA           = edpgi_self,
        Cognition    = edpgi_gencog_euro,
        BMI          = edpgi_bmi_euro,
        Alcoholism   = edpgi_alcohol_euro,
        Alzheimers   = edpgi_alzh_euro,
        Wellbeing    = edpgi_wellbeing_euro,
        Bipolar      = edpgi_bipolar_euro,
        ADHD         = edpgi_adhd_euro,
        Longevity    = edpgi_longevity_euro,
        Antisocial   = edpgi_antisocial_euro,
        Depression   = edpgi_mdepressived_euro,
        Anxiety      = edpgi_anxietyfactor_euro) %>%
    cor()
# Plot the matrix, only the lower left.
edpgi_correlation.matrix[
    upper.tri(edpgi_correlation.matrix, diag = FALSE)] <- NA
edpgi_correlation.plot <- edpgi_correlation.matrix %>%
    ggcorrplot::ggcorrplot(hc.order = FALSE,
        type = "full")
# Save this file
ggsave(file.path(figures.folder, "edpgi-correlation.png"),
    plot = edpgi_correlation.plot,
    units = "cm", width = 1.5 * fig.width, height = 1.25 * fig.height)
