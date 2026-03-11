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
library(fixest)
# Define number of digits in tables and graphs
digits.no <- 3
# Size for figures
fig.height <- 10
fig.width <- 1.25 * fig.height
presentation.width <- (3 / 2) * fig.width
presentation.height <- presentation.width
# List of 3 default colours.
colour.list <- c(
    "#1f77b4", # Blue
    "#2ca02c", # Green
    "#DAB1DA") #"#d62728") # Red
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

analysis.data %>% pull(edpgi_all_imputed_self) %>% summary() %>% print()
analysis.data %>% pull(edpgi_all_imputed_self) %>% sd() %>% print()
analysis.data %>% pull(edpgi_all_imputed_random) %>% summary() %>% print()
analysis.data %>% pull(edpgi_all_imputed_random) %>% var() %>% print()

# Get the demographic possible confounders.
demographic.data <- analysis.data %>%
    select(
        edpgi_all_imputed_self,
        edpgi_all_imputed_random,
        sex_male,
        recruitedage,
        sibling_count,
        urban,
        adhd_pgi,
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
    "t2diabetes_pgi"    = "Other PGI: Diabetes, Type-2",
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
        name = "Correlation Estimate, and 95% Confidence Interval",
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
        name = "Correlation Estimate",
        breaks = seq(-1, 1, by = 0.05)) +
    ggtitle(TeX("Demographic Information")) +
    theme(plot.title = element_text(size = rel(1), hjust = 0),
        plot.title.position = "plot",
        plot.margin = unit(c(0.5, 3, 0, 0), "mm"),
        legend.position = "bottom", #c(0.25, 0.95),
        axis.text.y = element_text(hjust = 0))

# Save this plot
ggsave(file.path(figures.folder, "demographic-correlates.png"),
    plot = combined.plot,
    dpi = 300, units = "cm",
    width = fig.width * 1.25, height = fig.height* 1.25)


################################################################################
## Validation of assumptions.

# Show the mean Ed PGI random component across the parent dist
analysis.data$parental_quantile <-
    ecdf(analysis.data$edpgi_all_imputed_parental)(
        analysis.data$edpgi_all_imputed_parental)
analysis.data$parental_quantile <- factor(
    round(9 * analysis.data$parental_quantile) + 1)
# Compare this to Houmark+ (2024) Figure 1.
analysis.data %>%
    group_by(parental_quantile) %>%
    summarise(
        edpgi_random_mean   = mean(edpgi_exclude_imputed_random),
        edpgi_random_sd     = sd(edpgi_exclude_imputed_random),
        edpgi_self_mean     = mean(edpgi_all_imputed_self),
        edpgi_self_sd       = sd(edpgi_all_imputed_self),
        edpgi_parental_mean = mean(edpgi_all_imputed_parental),
        edpgi_parental_sd   = sd(edpgi_all_imputed_parental)) %>%
    print()
mean(analysis.data$edpgi_exclude_imputed_random
    > analysis.data$edpgi_all_imputed_parental)

# Plot the distribution of self and random component.
library(ggridges)
# Ed PGI
parental_self.plot <- analysis.data %>%
    ggplot(aes(x = edpgi_all_imputed_self,
        y = parental_quantile)) +
    geom_density_ridges_gradient(
        scale = 3, rel_min_height = 0.01,
        colour = "black", fill = colour.list[2]) +
    geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
    theme_bw() +
    scale_x_continuous(expand = c(0, 0),
        name = "Ed PGI",
        breaks = seq(-10, 10, by = 1),
        limits = c(-4, 4)) +
    scale_y_discrete(expand = c(0, 0),
        name = "") +
    ggtitle("Parental Quantile, 1st through 10th") +
    theme(plot.title = element_text(size = rel(1), hjust = 0),
        plot.title.position = "plot",
        plot.margin = unit(c(0.5, 3, 0, 0), "mm"),
        panel.spacing = unit(0.1, "lines"),
        legend.position = "none")
# Save the plot.
ggsave(file.path(figures.folder, "edpgi-self-dist.png"),
    plot = parental_self.plot,
    units = "cm", width = fig.width, height = fig.height)

# Random component.
parental_random.plot <- analysis.data %>%
    ggplot(aes(x = edpgi_all_imputed_random,
        #x = (edpgi_all_imputed_self - edpgi_all_imputed_parental),
        y = parental_quantile)) +
    geom_density_ridges_gradient(
        scale = 3, rel_min_height = 0.01,
        colour = "black", fill = colour.list[3]) +
    geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
    theme_bw() +
    scale_x_continuous(expand = c(0, 0),
        name = "Ed PGI",
        breaks = seq(-10, 10, by = 1),
        limits = c(-4, 4)) +
    scale_y_discrete(expand = c(0, 0),
        name = "") +
    ggtitle("Parental Quantile, 1st through 10th") +
    theme(plot.title = element_text(size = rel(1), hjust = 0),
        plot.title.position = "plot",
        plot.margin = unit(c(0.5, 3, 0, 0), "mm"),
        panel.spacing = unit(0.1, "lines"),
        legend.position = "none")
# Save the plot.
ggsave(file.path(figures.folder, "edpgi-random-dist.png"),
    plot = parental_random.plot,
    units = "cm", width = fig.width, height = fig.height)


################################################################################
# Ensure this does not vary across the parental distribution.
analysis.data %>%
    lm(edpgi_all_imputed_self ~ 1 + edpgi_exclude_imputed_random * parental_quantile, data = .) %>%
    summary() %>%
    print()
#  TEST: does this association hold true among people for whom we observe one or both parents? -> Yes.
analysis.data %>%
    filter(father_present + mother_present > 0) %>%
    lm(edpgi_exclude_imputed_random ~ 1 + 
        + father_present : edpgi_all_imputed_paternal
        + mother_present : edpgi_all_imputed_maternal, data = .) %>%
    summary() %>%
    print()
# Compare to the Ed PGI in raw form.
analysis.data %>%
    filter(father_present + mother_present > 0) %>%
    lm(edpgi_all_imputed_self ~ 1 +
        + father_present : edpgi_all_imputed_paternal
        + mother_present : edpgi_all_imputed_maternal, data = .) %>%
    summary() %>%
    print()
