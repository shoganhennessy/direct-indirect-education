#!/usr/bin/R
## Senan Hogan-Hennessy, 24 March 2025
## UKB data -> First-stage analysis of PGI causal design (done with Ed PGI).
set.seed(47)
print(Sys.time())

# Functions for data manipulation and visualisation
library(tidyverse)
# Functions for tables into TeX
library(xtable)
# The standard, linear, IV estimator package.
library(ivreg)
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
    filter(analysis_sample == 1)
# Show the construction
print(analysis.data)
print(names(analysis.data))


################################################################################
## Testing out the first-stage, with random component of Ed PGI

# Correlation between random component and the realised PGI
analysis.data %>%
    lm(edpgi_self ~ 1 + edpgi_random, data = .) %>%
    summary() %>%
    print()
# Compare to the true DGP: PGI = E[ PGI | parents] + random component
analysis.data %>%
    lm(edpgi_self ~ 1 + edpgi_parents + edpgi_random, data = .) %>%
    summary() %>%
    print()
# Note coef on parents < 1 because of measurement error.
#TODO: use independent weightes for Ed PGI, obviously related IV robustness.


################################################################################
## Quasi second-stage, Ed PGI -> Ed Years

# OLS (not causal)
analysis.data %>%
    lm(edyears ~ 1 + edpgi_self, data = .) %>%
    summary() %>%
    print()

# Show the correlation between the random value, and Ed years
analysis.data %>%
    lm(edyears ~ 1 + edpgi_random, data = .) %>%
    summary() %>%
    print()

# Show the correlation between the self + parents value, and Ed years
analysis.data %>%
    lm(edyears ~ 1 + edpgi_self + edpgi_parents, data = .) %>%
    summary() %>%
    print()

# Try it with the random part as an instrument.
analysis.data %>%
    ivreg::ivreg(edyears ~ 1 + edpgi_self | 1 + edpgi_random,
        data = .) %>%
    summary() %>%
    print()

# Show the mean Ed PGI random component across the parent dist
analysis.data$parental_quantile <-
    ecdf(analysis.data$edpgi_parents)(analysis.data$edpgi_parents)
analysis.data$parental_quantile <- factor(
    round(10 * analysis.data$parental_quantile))
# Compare this to Houmark+ (2024) Figure 1.
analysis.data %>%
    group_by(parental_quantile) %>%
    summarise(
        edpgi_random_mean   = mean(edpgi_random),
        edpgi_random_sd     = sd(edpgi_random),
        edpgi_self_mean     = mean(edpgi_self),
        edpgi_self_sd       = sd(edpgi_self),
        edpgi_parental_mean = mean(edpgi_parents),
        edpgi_parental_sd   = sd(edpgi_parents)) %>%
    View()
mean(analysis.data$edpgi_self > analysis.data$edpgi_parents)
# Ensure this does not vary across the parental distribution.
analysis.data %>%
    lm(edpgi_self ~ 1 + edpgi_random * parental_quantile, data = .) %>%
    summary() %>%
    print()
#  TEST: does this association hold true among people for whom we observe one or both parents? -> Yes.
analysis.data %>%
    filter(father_present + mother_present > 0) %>%
    lm(edpgi_random ~ 1 + 
        sex_male * (edpgi_father + edpgi_mother), data = .) %>%
    summary() %>%
    print()
# Compare to the Ed PGI in raw form.
analysis.data %>%
    filter(father_present + mother_present > 0) %>%
    lm(edpgi_self ~ 1 +
        sex_male * (edpgi_father + edpgi_mother), data = .) %>%
    summary() %>%
    print()
