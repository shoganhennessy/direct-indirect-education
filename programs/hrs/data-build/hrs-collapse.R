#!/usr/bin/R
## Senan Hogan-Hennessy, 25 September 2023
## Collapse the pre-built HRS panel data set.
set.seed(47)
print(Sys.time())

# Functions for data manipulation and visualisation
library(tidyverse)
# Define number of digits in tables and graphs
digits.no <- 3
# Define folder path for the general data.
data.folder <- file.path("..", "..", "data", "hrs-public", "intermediate-files")
output.folder <- file.path("..", "..", "data", "hrs-public")


################################################################################
## Load HRS panel data

# Load the pre-cleaned HRS panel data.
hrs_panel.data <- data.folder %>%
    file.path("hrs-panel-cleaned.csv") %>%
    read_csv()

# Load the gene score data.
hrs_gene.data <- data.folder %>%
    file.path("hrs-gene-data.csv") %>%
    read_csv()

# Load the childhood SES data.
hrs_childhood.data <- data.folder %>%
    file.path("hrs-childhood-data.csv") %>%
    read_csv()


################################################################################
## Clean the panel data.

# Restrict to the relevant subsample.
hrs_panel.data <- hrs_panel.data %>%
    filter(
        # Those currently in the labour force, full-time working.
        lbrf == 1,
        # Those in final-working years, 50 <= age <= 65
        50 <= agey & agey <= 65
    ) %>%
    # Measure income in real (inflation adjusted) terms.
    mutate(
        indiv_earnings_real = iearn / cpi,
        indiv_income_real = itot / cpi)

# Make missing values from outlier incomes, income in (10, 300) x 10^3
hrs_panel.data <- hrs_panel.data %>%
    mutate(indiv_income_real = ifelse(
        10000 < indiv_income_real & indiv_income_real < 300000,
            indiv_income_real, NA)) %>%
    mutate(indiv_earnings_real = ifelse(
        10000 < indiv_earnings_real & indiv_earnings_real < 300000,
            indiv_earnings_real, NA))

# Names in the original data.
print(names(hrs_panel.data))

## Collapse the pre-cleaned HRS panel data, to get measure of LR income, etc.
# First, the demographic info
hrs_demographic.data <- hrs_panel.data %>%
    group_by(rahhidpn) %>%
    summarise(
        ## Demographic info:
        # Gender
        gender_female = max(as.numeric(ragender == 2), na.rm = TRUE),
        # Birth year rabyear
        indiv_birthyear = max(rabyear, na.rm = TRUE),
        # Race of individual
        race_white = min(as.numeric(raracem == 1), na.rm = TRUE),
        race_black = min(as.numeric(raracem == 2), na.rm = TRUE),
        race_other = min(as.numeric(raracem == 3), na.rm = TRUE),
        # maximum reported education years.
        indiv_edyears = max(raedyrs, na.rm = TRUE),
        # maximum reported degree obtained.
        degree_reported = max(raedegrm, na.rm = TRUE),
        # Parents' years of education.
        mother_edyears = max(rameduc, na.rm = TRUE),
        father_edyears = max(rafeduc, na.rm = TRUE),
        # Size of household, max ever reported.
        hhres = max(hhres, na.rm = TRUE),
        # Number of children, max ever reported.
        child = max(child, na.rm = TRUE),
        # Census division, max ever reported.
        cendiv = max(cendiv, na.rm = TRUE),
        # Physical BMI, max ever reported.
        pmbmi = max(pmbmi, na.rm = TRUE),
        # Physical height, max ever reported.
        pmhght = max(pmhght, na.rm = TRUE),
        # Physical weight, max ever reported.
        pmwght = max(pmwght, na.rm = TRUE),
        # Cognition score, best score measured
        cogtot = max(cogtot, na.rm = TRUE)) %>%
    ungroup()

# Second, info on income
hrs_income.data <- hrs_panel.data %>%
    filter(!is.na(indiv_earnings_real)) %>%
    group_by(rahhidpn) %>%
    summarise(
        ## Maximum reported income.
        # Number of years in the survey
        observed_count = n(),
        # maximum reported annual individual income.
        indiv_income_real = max(indiv_income_real, na.rm = TRUE),
        # maximum reported annual individual labour earnings
        indiv_earnings_real = max(indiv_earnings_real, na.rm = TRUE),
        # age at max income
        indiv_agey = agey[
            which(indiv_earnings_real  == max(indiv_earnings_real , na.rm = TRUE))],
        # year at max income
        survey_year = survey_year[
            which(indiv_earnings_real  == max(indiv_earnings_real , na.rm = TRUE))],
        # Survey weights, in year of max income
        wtresp = wtresp[
            which(indiv_earnings_real  == max(indiv_earnings_real , na.rm = TRUE))]) %>%
    ungroup()

# Second, info on gene scores
hrs_gene.data <- hrs_gene.data %>%
    filter(!is.na(genescore_educ_euro))

# Third info on childhood SES
hrs_childhood.data <- hrs_childhood.data %>%
    # Remove doubled up variables
    select(-rahhid, -father_edyears_supp, -mother_edyears_supp)

## Put all the three back together
hrs_collapsed.data <- hrs_demographic.data %>%
    left_join(hrs_income.data, by = "rahhidpn") %>%
    left_join(hrs_childhood.data, by = "rahhidpn") %>%
    left_join(hrs_gene.data, by = "rahhidpn")

# Save memory by removing the intermediary data files.
rm(list = c(
    "hrs_panel.data",
    "hrs_demographic.data",
    "hrs_income.data",
    "hrs_childhood.data",
    "hrs_gene.data"))
gc()

# Replace all Inf values with missings
# (these are generated by the summarise function on missing values.)
hrs_collapsed.data <- hrs_collapsed.data %>%
    mutate_if(is.double, list(~ na_if(., Inf))) %>%
    mutate_if(is.double, list(~ na_if(., -Inf)))

# Filter to non-missing in the Z, D, Y variables.
hrs_collapsed.data <- hrs_collapsed.data %>%
    filter(
        !is.na(genescore_educ_euro),
        !(is.na(degree_reported) & is.na(indiv_edyears)),
        !is.na(indiv_income_real))

# Restrict European ancestory, as gene score not reliable out of sample.
hrs_collapsed.data <- hrs_collapsed.data %>%
    filter(race_white == 1)

# Keep controls for parents' education level.
hrs_collapsed.data <- hrs_collapsed.data %>%
    mutate(parent_edyears = pmax(mother_edyears, father_edyears, na.rm = TRUE))


################################################################################
## Save the collapsed file.

# Save the collapsed panel data.
hrs_collapsed.data %>%
    write_csv(file.path(output.folder, "hrs-panel-collapsed.csv"))
