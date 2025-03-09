#!/usr/bin/R
## Senan Hogan-Hennessy, 2 July 2024
## Script to ingest raw HRS/RAND panel data, and put into a co
print(Sys.time())
set.seed(47)
# Packages:
# functions for data manipulation and visualisation
library(tidyverse)
# Package to load Stata files
library(haven)
# Define folder path for the raw data.
data.folder <- file.path("..", "..", "data")
output.folder <- file.path(data.folder, "hrs-public", "intermediate-files")


################################################################################
## Load raw data

#! Note: Can't run this part on a local machine -> run on the CCSS server.
# Load the combined RAND panel of HRS correspondents.
long.data <- data.folder %>%
    file.path("hrs-public",
        "randhrs1992_2020v2_STATA",
        "randhrs1992_2020v2.dta") %>%
    read_dta()


################################################################################
## Clean Longitudinal RAND HRS data.

# Select the relevant columns
clean_long.data <- long.data %>%
    # Create the person identifier.
    mutate(rahhid = str_sub(rahhidpn, 1, -4)) %>%
    # Restrict to the relevant columns
    select(
        ## Variables that DO NOT vary across the panel:
        # String identifiers, person + household
        rahhidpn, rahhid,
        # year of birth, Reference person
        rabyear,
        # year of death, Reference person
        radyear,
        # Gender, of reference person
        ragender,
        # Race, of reference person
        raracem,
        # Years of education, of reference person
        raedyrs,
        # Highest degree, of reference person
        raedegrm,
        # Mother's year of educ, of reference person
        rameduc,
        # Father's year of educ, of reference person
        rafeduc,
        # Whether a military veteran
        ravetrn,
        ## Variables that DO vary across the panel:
        # Person level weight
        ends_with("wtresp"),
        # Age at time of interview, reference person
        ends_with("agem_b"),
        ends_with("agey_b"),
        # Individual Earnings (both Ref + spouse, every wave of the survey).
        ends_with("iearn"),
        # Total household Earnings (household, every wave of the survey).
        ends_with("itot"),
        # Whether is part of a union.
        ends_with("union"),
        # Labour force status (both Ref + spouse, every wave of the survey).
        ends_with("lbrf"),
        # Whether retired (both Ref + spouse, every wave of the survey).
        ends_with("sayret"),
        # HH num of people
        ends_with("hhres"),
        # HH num of children
        ends_with("child"),
        # Census division
        ends_with("cendiv"),
        # Physical measures, BMI
        ends_with("pmbmi"),
        # Physical measures, height
        ends_with("pmhght"),
        # Physical measures, weight
        ends_with("pmwght"),
        # Reasoning scores, total score for cognitive status.
        ends_with("cogtot"),
        # Reasoning scores, number series score
        ends_with("nsscre"),
        # Reasoning scores, vocab series score
        ends_with("vocab"),
        # Reasoning scores, total word recall score
        ends_with("tr20"),
        # Reasoning scores, immediate word recall
        ends_with("imrc"),
        # Reasoning scores, delayed word recall
        ends_with("dlrc"),
        # Reasoning scores, backwards word counting
        ends_with("bwc20"),
        # Reasoning scores, able to name president
        ends_with("pres"))

# Designate the identifying variables
nonvarying.list <- clean_long.data %>%
    select(starts_with("ra")) %>%
    names()

# Convert to a panel (and out of a wide format).
hrs_panel.data <- clean_long.data %>%
    pivot_longer(
        # Designate the identifying variables
        cols = - nonvarying.list,
        # Naming patter: RWvar, where
        # R \in {h,r,s}, W = 1,...,15 HRS wave, var is string for variable name.
        names_to = c("person", "wave", "variable"),
        names_pattern = "([hrs])([0-9]+)([A-Za-z]+)",
        values_to = "value",
        values_drop_na = FALSE) %>%
    # Restrict to household reference persons (and not spouse)
    filter(person != "s") %>%
    # Recode household variables to the reference person
    mutate(person = ifelse(person == "h", "r", person)) %>%
    # Make a panel format (i.e., time varying variables as columns)
    pivot_wider(
        id_cols = c("person", "wave", nonvarying.list),
        names_from = variable, values_from = value)

# Get the year of the survey, based on wave
hrs_panel.data <- hrs_panel.data %>%
    mutate(survey_year =
        ifelse(is.na(wave), NA,
        ifelse(wave == 1,   1992,
        ifelse(wave == 2,   1994,
        ifelse(wave == 3,   1996,
        ifelse(wave == 4,   1998,
        ifelse(wave == 5,   2000,
        ifelse(wave == 6,   2002,
        ifelse(wave == 7,   2004,
        ifelse(wave == 8,   2006,
        ifelse(wave == 9,   2008,
        ifelse(wave == 10,  2010,
        ifelse(wave == 11,  2012,
        ifelse(wave == 12,  2014,
        ifelse(wave == 13,  2016,
        ifelse(wave == 14,  2018,
        ifelse(wave == 15,  2020, NA)))))))))))))))))


################################################################################
## Save the resulting data file.

# Save the gene score data file.
hrs_panel.data %>%
    write_csv(file.path(output.folder, "hrs-panel-data.csv"))
