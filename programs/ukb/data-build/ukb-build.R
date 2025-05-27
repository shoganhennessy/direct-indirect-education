#!/usr/bin/R
## Senan Hogan-Hennessy, 11 March 2025.
## Script to ingest raw UKB data, and output familiar tabular format.
print(Sys.time())
set.seed(47)
## Packages:
# functions for data manipulation and visualisation
library(tidyverse)
# Functions for fast data manipulation (adapting to data.table code from Kweon+)
library(data.table)

# Define folder paths (1) input data (2) clean data.
data.folder <- file.path("..", "..", "..", "data", "ukb-restricted")
input.folder <- file.path(data.folder, "input")
output.folder <- file.path(data.folder, "cleaned")


################################################################################
## Load raw data

# Load the phenotype data (raw format from UKB RAP servers).
ukb_pheno.data <- input.folder %>%
    file.path("phenotype-extract.csv") %>%
    read_csv()

# SHow how these data look.
print(ukb_pheno.data)
print(names(ukb_pheno.data))

# Load the relative connections data
ukb_relatives.data <- input.folder %>%
    file.path("kinship-adjusted.dat") %>%
    read_tsv()

# Load the Ed PGI (Okbay+ 2022) data, raw format from UKB RAP servers.
ukb_edpgi.data <- input.folder %>%
    file.path("ed-pgi-score.tsv") %>%
    read_tsv()

# Load the imputed Ed PGI (Young+ 2022), raw format from UKB RAP servers.
ukb_imputed.data <- input.folder %>%
    file.path("imputed-ed-pgi.pgs.txt") %>%
    read_table()

# Load the occupation-coded income data (thanks to Kweon Koellinger+ 2025).
load(file.path(input.folder, "earnings-imputed", "data_input.Rdata"))

# CPI-H (including housing costs) for the UK (2015 == 100).
# https://www.ons.gov.uk/economy/inflationandpriceindices/datasets/consumerpriceindices
cpi.data <- input.folder %>%
    file.path("cpih-uk.csv") %>%
    read_csv()


################################################################################
## Clean the multiple instance columns, and standardise ed quals variables.

# Code for making the instance files homogenised.
ukb_pheno.data <- ukb_pheno.data %>%
    mutate(
        # SOC job code
        jobcode_soc =
            ifelse(!is.na(participant.p20277_i0), participant.p20277_i0,
            ifelse(!is.na(participant.p20277_i1), participant.p20277_i1,
            ifelse(!is.na(participant.p20277_i2), participant.p20277_i2,
            ifelse(!is.na(participant.p20277_i3), participant.p20277_i3, NA)))),
        # birth_n_coord
        birth_n_coord =
            ifelse(!is.na(participant.p129_i0) & participant.p129_i0 != -1, participant.p129_i0,
            ifelse(!is.na(participant.p129_i1) & participant.p129_i1 != -1, participant.p129_i1,
            ifelse(!is.na(participant.p129_i2) & participant.p129_i2 != -1, participant.p129_i2, NA))),
        # birth_e_coord
        birth_e_coord =
            ifelse(!is.na(participant.p130_i0), participant.p130_i0,
            ifelse(!is.na(participant.p130_i1), participant.p130_i1,
            ifelse(!is.na(participant.p130_i2), participant.p130_i2, NA))),
        # birth_country
        birth_country =
            ifelse(!is.na(participant.p1647_i0), participant.p1647_i0,
            ifelse(!is.na(participant.p1647_i1), participant.p1647_i1,
            ifelse(!is.na(participant.p1647_i2), participant.p1647_i2, NA))),
        # employment
        employment =
            ifelse(!is.na(participant.p6142_i0), participant.p6142_i0,
            ifelse(!is.na(participant.p6142_i1), participant.p6142_i1,
            ifelse(!is.na(participant.p6142_i2), participant.p6142_i2,
            ifelse(!is.na(participant.p6142_i3), participant.p6142_i3, NA)))),
        # hours_workweek
        hours_workweek =
            ifelse(!is.na(participant.p767_i0) & participant.p767_i0 > 0, participant.p767_i0,
            ifelse(!is.na(participant.p767_i1) & participant.p767_i0 > 0, participant.p767_i1,
            ifelse(!is.na(participant.p767_i2) & participant.p767_i0 > 0, participant.p767_i2,
            ifelse(!is.na(participant.p767_i3) & participant.p767_i0 > 0, participant.p767_i3,
                NA)))))

# Code edQuals -> Edyears Following ISCED
ukb_pheno.data <- ukb_pheno.data %>%
    mutate(edyears =
        ifelse(participant.edqual_highered == 1,     18, # Uni degree -> 18
        ifelse(participant.edqual_professional == 1, 15, # Professional post-secondary -> 15
        ifelse(participant.edqual_alevels == 1,      14, # A Levels -> 14
        ifelse(participant.edqual_vocational == 1,   14, # Vocational post-secondary -> 14 
        ifelse(participant.edqual_gcses == 1,        12, # GCSEs -> 12
        ifelse(participant.edqual_minimum == 1,      9,  # None of above -> lower secondary (i.e., legal minimum). 
        ifelse(participant.edqual_missing == 1,      NA, # Missing -> NA
            100)))))))) # End condition for everyone.
# Ensure reasonable distribution for everyone (and end condition not met.)
ukb_pheno.data %>% pull(edyears) %>% table(exclude = NULL) %>% print()
ukb_pheno.data %>% pull(jobcode_soc) %>% table(exclude = NULL) %>% print()


################################################################################
## Clean Pheno file.

# Clean the phenotype data.
cleaned_pheno.data <- ukb_pheno.data %>%
    # Restrict data types
    transmute(
        eid                    = as.integer(participant.eid),
        sex_male               = as.integer(participant.p31),
        visitdate              = date(participant.p53_i0),
        deathdate              = date(participant.p40000_i0),
        urban_cat              = as.integer(participant.p20118_i0),
        # Birth info.
        birth_n_coord          = as.integer(birth_n_coord),
        birth_e_coord          = as.integer(birth_e_coord),
        birth_country          = as.integer(birth_country),
        birthyear              = as.integer(participant.p34),
        birthmonth             = as.integer(participant.p52),
        # Work + income categories,
        householdincome_cat    = as.integer(participant.p738_i0),
        recruitedage           = as.integer(participant.p21022),
        jobcode_soc            = as.integer(jobcode_soc),
        employment             = as.integer(employment),
        hours_workweek         = as.integer(hours_workweek),
        # Education info.
        edqual_highered        = as.integer(participant.edqual_highered),
        edqual_professional    = as.integer(participant.edqual_professional),
        edqual_alevels         = as.integer(participant.edqual_alevels),
        edqual_vocational      = as.integer(participant.edqual_vocational),
        edqual_gcses           = as.integer(participant.edqual_gcses),
        edqual_minimum         = as.integer(participant.edqual_minimum),
        edqual_missing         = as.integer(participant.edqual_missing),
        edyears                = as.integer(edyears),
        # Genetic variables.
        genetic_race           = as.integer(participant.p22006),
        PCA                    = as.integer(participant.p22020),
        asthma_pgi             = as.integer(participant.p26210),  
        bipolar_pgi            = as.integer(participant.p26214),  
        bmi_pgi                = as.integer(participant.p26216),  
        height_pgi             = as.integer(participant.p26240),  
        schizophrenia_pgi      = as.integer(participant.p26275),  
        t2diabetes_pgi         = as.integer(participant.p26285)
        ) %>%
    # Clean the resulting columns
    mutate(
        # householdincome_cat,    p738_i0,    https://biobank.ctsu.ox.ac.uk/crystal/field.cgi?id=738
        householdincome_less18k   = ifelse(is.na(householdincome_cat), NA,
            as.integer(householdincome_cat == 1)),
        householdincome_18to31k   = ifelse(is.na(householdincome_cat), NA,
            as.integer(householdincome_cat == 2)),
        householdincome_31to52k   = ifelse(is.na(householdincome_cat), NA,
            as.integer(householdincome_cat == 3)),
        householdincome_52to100k  = ifelse(is.na(householdincome_cat), NA,
            as.integer(householdincome_cat == 4)),
        householdincome_above100k = ifelse(is.na(householdincome_cat), NA,
            as.integer(householdincome_cat == 5)),
        # Replace missing values with binary values.
        genetic_euroancest = ifelse(is.na(genetic_race), 0, genetic_race),
        inPCA = ifelse(is.na(PCA), 0, PCA),
        # death date variables
        deathyear = as.integer(format(deathdate, "%Y")),
        visityear = as.integer(format(visitdate, "%Y")),
        hasdied = as.integer(!is.na(deathdate)))


################################################################################
## Connect UKB data with collected distance to nearest uni data.

# Read the birth counties data.
#TODO might not use this one.
birth_counties.data <- read_csv(
    file.path(input.folder, "birth-locations", "main_phase_locations.csv"))

# Read the higher ed locations file.
higher_loc.data <- read_csv(
    file.path(input.folder, "..", "..", "uk-highered", "highered-compiled.csv"))

#TODO: use the counties base map to assign counties to 
#TODO: birth_n_coord, birth_e_coord

#TODO: write merging code that gives the closest uni
#TODO: (if currently founded in year aged 17).
#TODO: Then a binary for whether there is a uni in your county, and distance to it.


################################################################################
## Clean PGI data files.

# Restrict Ed PGI to identifiers, and the score
cleaned_edpgi.data <- ukb_edpgi.data %>%
    transmute(
        eid = IID,
        edpgi_all_raw = SCORE1_AVG)

# Restrict imputed variables identifiers, and the score
cleaned_imputed.data <- ukb_imputed.data %>%
    transmute(
        eid                        = IID,
        famid                      = FID,
        edpgi_imputed_self_raw     = proband,
        edpgi_imputed_sibling_raw  = sibling,
        edpgi_imputed_paternal_raw = paternal,
        edpgi_imputed_maternal_raw = maternal)

#TODO: get a second PGI measure for Ed PGI, to use ORIV.


################################################################################
## Calculate count of family members', and whether their parents are in data.

# Note: relative connections are by-directional, so need to append to self. 
count_relatives.data <- ukb_relatives.data %>%
    # Note: relative data are one-directional,
    #       -> need to append to self to count second column, too. 
    rename(ID1 = ID2, ID2 = ID1) %>%
    bind_rows(ukb_relatives.data) %>%
    mutate(
        eid = ID1,
        sibling = as.integer(InfType == "FS"),
        parent = as.integer(InfType == "PO")) %>%
    # Count parent + sibling connections in these data.
    group_by(eid) %>%
    summarise(
        sibling_count = sum(sibling, na.rm = TRUE),
        parent_count = sum(parent, na.rm = TRUE)) %>%
    ungroup()

# Column to show whether mother + father are in the data.
relatives_present.data <- ukb_imputed.data %>%
    transmute(
        eid = IID,
        father_present  = as.integer(!is.na(FATHER_ID)),
        mother_present  = as.integer(!is.na(MOTHER_ID)),
        sibling_present = as.integer(!is.na(sibling)))

# Show how many parents there are in the total data.
print(sum(relatives_present.data$father_present))
print(sum(relatives_present.data$mother_present))
print(sum(relatives_present.data$father_present) +
    sum(relatives_present.data$mother_present))


################################################################################
## Calculate OCC coded wages, thanks to Kweon Koellinger provided code.

# Get the data file to merge with OCC coded wage data.
ukb_soc.data <- cleaned_pheno.data %>%
    filter(!is.na(jobcode_soc)) %>%
    # Get the variables needed for income-impute.R
    transmute(eid = eid,
        SOC  = as.integer(jobcode_soc),
        male = as.integer(sex_male),
        age  = as.integer(recruitedage),
        year = as.integer(visityear))

# Use the Kweon Koellinger script for SOC occ hourly wages.
source("income-impute.R")
ukb_soc.data <- GET_INCOME(ukb_soc.data)

# Update to 2024 (income-impute gives 2015 GBP, so update to year specified.)
year.base <- 2015
cpi.base <- cpi.data %>% filter(year == year.base) %>% pull(cpi2015)
year.update <- 2024
cpi.update <- cpi.data %>% filter(year == year.update) %>% pull(cpi2015)
cpi.factor <- cpi.update / cpi.base
# Adjust the columns.
ukb_soc.data <- ukb_soc.data %>%
    transmute(
        eid = eid,
        soc_mean_hourly       = cpi.factor * soc_mean_hourly)
        # soc_median_hourly     = cpi.factor * soc_median_hourly,
        # soc_mean_hourly_all   = cpi.factor * soc_mean_hourly_all,
        # soc_median_hourly_all = cpi.factor * soc_median_hourly_all,
        # soc_2d_hourly         = exp(log_y_hourly))


################################################################################
## Merge data files.

# Merge the Phenotype data with Ed PGI.
cleaned_pheno.data <- cleaned_pheno.data %>%
    left_join(cleaned_edpgi.data, by = "eid")

# Merge the Phenotype data with parental imputed Ed PGI.
cleaned_pheno.data <- cleaned_pheno.data %>%
    left_join(cleaned_imputed.data, by = "eid")

# Get counts of relatives onto phenotype data.
cleaned_pheno.data <- cleaned_pheno.data %>%
    left_join(count_relatives.data,   by = "eid") %>%
    left_join(relatives_present.data, by = "eid") %>%
    # Missings are zeros.
    mutate(
        sibling_count   = ifelse(is.na(sibling_count),   0, sibling_count),
        parent_count    = ifelse(is.na(parent_count),    0, parent_count),
        father_present  = ifelse(is.na(father_present),  0, father_present),
        mother_present  = ifelse(is.na(mother_present),  0, mother_present),
        sibling_present = ifelse(is.na(sibling_present), 0, sibling_present))

# Wage data onto phenotype data.
cleaned_pheno.data <- cleaned_pheno.data %>%
    left_join(ukb_soc.data, by = "eid")


################################################################################
## Clean the Ed PGI data fields, within entire sample.

# Standardise the Ed PGI.
cleaned_pheno.data <- cleaned_pheno.data %>%
    # PGI reliable for Euro ancestry, and reliable PCA data.
    mutate(
        edpgi_all_norm =
            ifelse(genetic_euroancest == 0, NA, edpgi_all_raw),
        edpgi_imputed_self_norm =
            ifelse(genetic_euroancest == 0, NA, edpgi_imputed_self_raw),
        edpgi_imputed_sibling_norm =
            ifelse(genetic_euroancest == 0, NA, edpgi_imputed_sibling_raw),
        edpgi_imputed_paternal_norm =
            ifelse(genetic_euroancest == 0, NA, edpgi_imputed_paternal_raw),
        edpgi_imputed_maternal_norm =
            ifelse(genetic_euroancest == 0, NA, edpgi_imputed_maternal_raw)
    ) %>%
    # Normalise to have 0 mean, 1 SD.
    # The parent/sibling values are centred around proband mean, sd
    mutate(
        edpgi_all_norm = (edpgi_all_norm -
            mean(edpgi_all_norm, na.rm = TRUE)) / sd(edpgi_all_norm, na.rm = TRUE),
        edpgi_imputed_self_norm = (edpgi_imputed_self_norm -
            mean(edpgi_imputed_self_norm, na.rm = TRUE)) / sd(edpgi_imputed_self_norm, na.rm = TRUE),
        edpgi_imputed_sibling_norm = (edpgi_imputed_sibling_norm -
            mean(edpgi_imputed_self_norm, na.rm = TRUE)) / sd(edpgi_imputed_self_norm, na.rm = TRUE),
        edpgi_imputed_paternal_norm = (edpgi_imputed_paternal_norm -
            mean(edpgi_imputed_self_norm, na.rm = TRUE)) / sd(edpgi_imputed_self_norm, na.rm = TRUE),
        edpgi_imputed_maternal_norm = (edpgi_imputed_maternal_norm -
            mean(edpgi_imputed_self_norm, na.rm = TRUE)) / sd(edpgi_imputed_self_norm, na.rm = TRUE))

#! Test: did the rescaling work?
cleaned_pheno.data %>%
    summarise(
        self_mean    = mean(edpgi_imputed_self_norm, na.rm = TRUE),
        self_sd      = sd(edpgi_imputed_self_norm, na.rm = TRUE),
        sibling_mean = mean(edpgi_imputed_sibling_norm, na.rm = TRUE),
        sibling_sd   = sd(edpgi_imputed_sibling_norm, na.rm = TRUE),
        father_mean  = mean(edpgi_imputed_paternal_norm, na.rm = TRUE),
        father_sd    = sd(edpgi_imputed_paternal_norm, na.rm = TRUE),
        mother_mean  = mean(edpgi_imputed_maternal_norm, na.rm = TRUE),
        mother_sd    = sd(edpgi_imputed_maternal_norm, na.rm = TRUE)) %>%
    print()

# Ensure outliers in PGIs do not persist.
hist(cleaned_pheno.data$edpgi_all_raw)
hist(cleaned_pheno.data$edpgi_all_norm)
# Similarly for imputed version.
hist(cleaned_pheno.data$edpgi_imputed_self_raw)
hist(cleaned_pheno.data$edpgi_imputed_self_norm)

# Show observation counts by both methods.
cleaned_pheno.data %>% filter(!is.na(edpgi_imputed_self_norm)) %>% print()
cleaned_pheno.data %>% filter(!is.na(edpgi_all_norm)) %>% print()

# Validate the correlation between two different packages for computing Ed PGI.
print(cor(cleaned_pheno.data$edpgi_all_norm,
    cleaned_pheno.data$edpgi_imputed_self_norm,
    use = "pairwise.complete.obs"))


################################################################################
## Declare base data file.

# Rename and restrict to necessary columns.
final_pheno.data <- cleaned_pheno.data %>%
    mutate(
        # Ed PGI among entire sample
        edpgi_all     = edpgi_all_norm,
        # Ed PGI among sample with siblings -> main PGI in my analysis
        edpgi_self    = edpgi_imputed_self_norm,
        edpgi_father  = edpgi_imputed_paternal_norm,
        edpgi_mother  = edpgi_imputed_maternal_norm,
        # Annual wage by hours worked
        soc_mean_annual = 52.14 * soc_mean_hourly * hours_workweek) %>%
    # Mean parental Ed PGI, scaled by 
    rowwise() %>%
    mutate(edpgi_parents = mean(
        c(edpgi_imputed_paternal_norm, edpgi_imputed_maternal_norm),
            na.rm = TRUE)) %>%
    ungroup() %>%
    tibble() %>%
    # Mark the analysis sample, those with imputed Ed PGI + Ed + SOC data.
    mutate(analysis_sample = as.integer(
        !is.na(edpgi_parents) & !is.na(edyears) & !is.na(soc_mean_annual))) %>%
    # Select the relevant columns.
    select(eid, famid,
        # Demographic variables.
        sex_male,
        visityear,
        recruitedage,
        genetic_euroancest,
        urban_cat,
        # Birth variables.
        birth_n_coord,
        birth_e_coord,
        birth_country,
        birthyear,
        birthmonth,
        # Genetic variables.
        edpgi_all,
        edpgi_self,
        edpgi_parents,
        edpgi_father,
        edpgi_mother,
        asthma_pgi,
        bipolar_pgi,
        bmi_pgi,
        height_pgi,
        schizophrenia_pgi,
        t2diabetes_pgi,
        # Education variables
        edqual_highered,
        edqual_professional,
        edqual_alevels,
        edqual_vocational,
        edqual_gcses,
        edqual_minimum,
        edqual_missing,
        edyears,
        # Income variables
        householdincome_cat,
        recruitedage,
        jobcode_soc,
        #TODO employment -> fix in the python extract because lists.
        hours_workweek,
        soc_mean_hourly,
        soc_mean_annual,
        starts_with("householdincome"),
        # Health & Death variables
        deathyear,
        # Designators
        sibling_count,
        parent_count,
        father_present,
        mother_present,
        sibling_present,
        analysis_sample)

# Show how large the relevant data fields are.
final_pheno.data %>% NROW()
final_pheno.data %>% filter(analysis_sample == 1) %>% NROW()

# Similarly, for with SOC occ wage data.
final_pheno.data %>% filter(!is.na(soc_mean_hourly)) %>% NROW()
final_pheno.data %>% filter(!is.na(edpgi_parents) & !is.na(edyears), !is.na(soc_mean_hourly)) %>% NROW()


################################################################################
## Calculate the random component of Ed PGI.

# Get the sample for whom we are getting mean expected PGI.
analysis.data <- final_pheno.data %>%
    filter(!is.na(edpgi_parents))

# (1) Calculate E[ Ed PGI | parents]
expected_pgi.reg <- analysis.data %>%
    lm(edpgi_self ~ poly(edpgi_father, 3) * poly(edpgi_mother, 3
        ) * sex_male * father_present * mother_present,
        data = .)
# (2) Calculate random component = Ed PGI - \hat E[ Ed PGI | parents]
analysis.data$edpgi_random <-
    analysis.data$edpgi_self - predict(expected_pgi.reg, analysis.data)

# Show summary statistics of the random component.
print(mean(analysis.data$edpgi_random))
print(sd(analysis.data$edpgi_random))
hist(analysis.data$edpgi_random)

# Put Ed PGi random component back onto the phenotype data file.
final_pheno.data <- analysis.data %>%
    select(eid, edpgi_random) %>%
    right_join(final_pheno.data, by = "eid")


################################################################################
## Save the resulting data file.

# Save the file.
final_pheno.data %>%
    write_csv(file.path(output.folder, "ukb-cleaned-pheno.csv"))
