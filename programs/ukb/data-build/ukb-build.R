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
# Geospatial filetypes.
library(sf)

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

# Load the occupation-coded income data (thanks to Kweon Koellinger+ 2025).
load(file.path(input.folder, "earnings-imputed", "data_input.Rdata"))

# CPI-H (including housing costs) for the UK (2015 == 100).
# https://www.ons.gov.uk/economy/inflationandpriceindices/datasets/consumerpriceindices
cpi.data <- input.folder %>%
    file.path("cpih-uk.csv") %>%
    read_csv()

## Load each of the PGI data files.

# Load the relative connections data
ukb_relatives.data <- input.folder %>%
    file.path("kinship-adjusted.dat") %>%
    read_tsv()

## Okbay+ (2022) Ed PGI, among the 23+Me sample.
# Load the Ed PGI (Okbay+ 2022) data, raw format from UKB RAP servers.
ukb_raw_edpgi_all.data <- input.folder %>%
    file.path("pgi", "pgi-okbay-all-2022", "raw-ed-pgi-okbay-all-2022.tsv") %>%
    read_tsv()

# Load the imputed Ed PGI (Young+ 2022), raw format from UKB RAP servers.
ukb_imputed_edpgi_all.data <- input.folder %>%
    file.path("pgi", "pgi-okbay-all-2022",
        "imputed-ed-pgi-okbay-all-2022.pgs.txt") %>%
    read_table()

## Okbay+ (2022) Ed PGI, excluding the 23+Me sample.
# Load the Ed PGI (Okbay+ 2022) data, raw format from UKB RAP servers.
ukb_raw_edpgi_exclude.data <- input.folder %>%
    file.path("pgi", "pgi-okbay-exclude-2022",
        "raw-ed-pgi-okbay-exclude-2022.tsv") %>%
    read_tsv()

# Load the imputed Ed PGI (Young+ 2022), raw format from UKB RAP servers.
ukb_imputed_edpgi_exclude.data <- input.folder %>%
    file.path("pgi", "pgi-okbay-exclude-2022",
        "imputed-ed-pgi-okbay-exclude-2022.pgs.txt") %>%
    read_table()

## ADHD PGI, hand computed (since not provided in UKB pheno data).
ukb_raw_adhdpgi.data <- input.folder %>%
    file.path("pgi", "ADHD", "raw-adhd-pgi.tsv") %>%
    read_tsv()


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

# Recode inconsistencies in edquals.
# GCSEs -> A-Levels -> Higher ed are sequentially dependent, so should be so.
ukb_pheno.data <- ukb_pheno.data %>%
    mutate(
        # if highered == 1, then A-Levels and GCSES == 1, by necessity
        participant.edqual_alevels = ifelse(participant.edqual_highered == 1,
            1, participant.edqual_alevels),
        participant.edqual_gcses = ifelse(participant.edqual_highered == 1,
            1, participant.edqual_gcses),
        # if A-Levels == 1, then GCSES == 1, by necessity
        participant.edqual_gcses = ifelse(participant.edqual_alevels == 1,
            1, participant.edqual_gcses))

#TODO: assign education years as the age they reported leaving full-time 
#TODO: education minus five for those with a vocational degree
#TODO: as their highest qualification.  

# Code edQuals -> Edyears Following ISCED
ukb_pheno.data <- ukb_pheno.data %>%
    mutate(edyears =
        ifelse(participant.edqual_highered == 1,     18, # Uni degree -> 18
        ifelse(participant.edqual_professional == 1, 11, # Professional post-secondary -> 15
        ifelse(participant.edqual_alevels == 1,      14, # A Levels -> 14
        ifelse(participant.edqual_vocational == 1,   11, # Vocational post-secondary -> 14 
        ifelse(participant.edqual_gcses == 1,        12, # GCSEs -> 12
        ifelse(participant.edqual_minimum == 1,      9,  # None of above -> lower secondary (i.e., legal minimum). 
        ifelse(participant.edqual_missing == 1,      NA, # Missing -> NA
            100)))))))) # End condition for everyone.
# Ensure reasonable distribution for everyone (and NA end condition not met.)
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
        #TODO: adjust python for lists.employment             = as.integer(employment),
        hours_workweek         = as.numeric(hours_workweek),
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
        genetic_race           = as.numeric(participant.p22006),
        PCA                    = as.numeric(participant.p22020),
        asthma_pgi             = as.numeric(participant.p26210),
        bipolar_pgi            = as.numeric(participant.p26214),
        bmi_pgi                = as.numeric(participant.p26216),
        height_pgi             = as.numeric(participant.p26240),
        schizophrenia_pgi      = as.numeric(participant.p26275),
        t2diabetes_pgi         = as.numeric(participant.p26285)
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
        householdincome_missing = ifelse(is.na(householdincome_cat) |
            householdincome_cat == 6 | householdincome_cat == 7, 1, 0),
        # Replace missing values with binary values.
        genetic_euroancest = ifelse(is.na(genetic_race), 0, genetic_race),
        inPCA = ifelse(is.na(PCA), 0, PCA),
        # death date variables
        deathyear = as.integer(format(deathdate, "%Y")),
        visityear = as.integer(format(visitdate, "%Y")),
        hasdied = as.integer(!is.na(deathdate)))


################################################################################
## Clean PGI data files.

# The raw PGIs, for (nearly) all 500,000 UKB people
# Restrict Ed PGI to identifiers, and the score
ukb_raw_edpgi_all.data <- ukb_raw_edpgi_all.data %>%
    transmute(
        eid = IID,
        edpgi_all_raw_self = SCORE1_AVG)

ukb_raw_edpgi_exclude.data <- ukb_raw_edpgi_exclude.data %>%
    transmute(
        eid = IID,
        edpgi_exclude_raw_self = SCORE1_AVG)

ukb_raw_adhdpgi.data <- ukb_raw_adhdpgi.data %>%
    transmute(
        eid = IID,
        adhd_pgi = SCORE1_AVG)

## The two imputed Ed PGIs, for the sibling subsample.
# Get the columns referring to own PGI, and family members' 
ukb_imputed_edpgi_all.data <- ukb_imputed_edpgi_all.data  %>%
    transmute(
        eid                        = IID,
        famid                      = FID,
        edpgi_all_imputed_self     = proband,
        #edpgi_all_imputed_sibling  = sibling,
        edpgi_all_imputed_paternal = paternal,
        edpgi_all_imputed_maternal = maternal)

ukb_imputed_edpgi_exclude.data <- ukb_imputed_edpgi_exclude.data %>%
    transmute(
        eid                            = IID,
        #famid                          = FID,
        edpgi_exclude_imputed_self     = proband,
        #edpgi_exclude_imputed_sibling  = sibling,
        edpgi_exclude_imputed_paternal = paternal,
        edpgi_exclude_imputed_maternal = maternal)


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
relatives_present.data <- input.folder %>%
    file.path("pgi", "pgi-okbay-all-2022",
        "imputed-ed-pgi-okbay-all-2022.pgs.txt") %>%
    read_table() %>%
    transmute(
        eid = IID,
        father_present  = as.integer(!is.na(FATHER_ID)),
        mother_present  = as.integer(!is.na(MOTHER_ID)),
        sibling_present = as.integer(!is.na(sibling)))

# Show how many parents there are in the total data.
print(sum(relatives_present.data$father_present))
print(sum(relatives_present.data$mother_present))
print(sum((relatives_present.data$father_present +
    relatives_present.data$mother_present) > 1))


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

# Merge the ADHD PGI.
cleaned_pheno.data <- cleaned_pheno.data %>%
    left_join(ukb_raw_adhdpgi.data, by = "eid")

# Merge the Ed PGIs (two different ones), both raw + imputed.
cleaned_pheno.data <- cleaned_pheno.data %>%
    left_join(ukb_raw_edpgi_all.data, by = "eid") %>%
    left_join(ukb_raw_edpgi_exclude.data, by = "eid") %>%
    left_join(ukb_imputed_edpgi_all.data, by = "eid") %>%
    left_join(ukb_imputed_edpgi_exclude.data, by = "eid")

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

# PGIs only reliable among the European ancestory (unfortunately).
genetic_euro_NA <- ifelse(cleaned_pheno.data$genetic_euroancest == 0, NA,
    cleaned_pheno.data$genetic_euroancest)
cleaned_pheno.data$asthma_pgi <- (genetic_euro_NA * cleaned_pheno.data$asthma_pgi)
cleaned_pheno.data$bipolar_pgi <- (genetic_euro_NA * cleaned_pheno.data$bipolar_pgi)                   
cleaned_pheno.data$bmi_pgi <- (genetic_euro_NA * cleaned_pheno.data$bmi_pgi)
cleaned_pheno.data$height_pgi <- (genetic_euro_NA * cleaned_pheno.data$height_pgi)                    
cleaned_pheno.data$schizophrenia_pgi <- (genetic_euro_NA * cleaned_pheno.data$schizophrenia_pgi)
cleaned_pheno.data$t2diabetes_pgi <- (genetic_euro_NA * cleaned_pheno.data$t2diabetes_pgi)                
cleaned_pheno.data$adhd_pgi <- (genetic_euro_NA * cleaned_pheno.data$adhd_pgi)
cleaned_pheno.data$edpgi_all_raw_self <- (genetic_euro_NA * cleaned_pheno.data$edpgi_all_raw_self)            
cleaned_pheno.data$edpgi_exclude_raw_self <- (genetic_euro_NA * cleaned_pheno.data$edpgi_exclude_raw_self)
cleaned_pheno.data$edpgi_all_imputed_self <- (genetic_euro_NA * cleaned_pheno.data$edpgi_all_imputed_self)
cleaned_pheno.data$edpgi_all_imputed_paternal <- (genetic_euro_NA * cleaned_pheno.data$edpgi_all_imputed_paternal)
cleaned_pheno.data$edpgi_all_imputed_maternal <- (genetic_euro_NA * cleaned_pheno.data$edpgi_all_imputed_maternal)    
cleaned_pheno.data$edpgi_exclude_imputed_self <- (genetic_euro_NA * cleaned_pheno.data$edpgi_exclude_imputed_self)
cleaned_pheno.data$edpgi_exclude_imputed_paternal <- (genetic_euro_NA * cleaned_pheno.data$edpgi_exclude_imputed_paternal)
cleaned_pheno.data$edpgi_exclude_imputed_maternal <- (genetic_euro_NA * cleaned_pheno.data$edpgi_exclude_imputed_maternal)

# Scale parental values by the children SD
cleaned_pheno.data$edpgi_all_imputed_paternal <- (cleaned_pheno.data$edpgi_all_imputed_paternal
    - mean(cleaned_pheno.data$edpgi_all_imputed_self, na.rm = TRUE)) / sd(cleaned_pheno.data$edpgi_all_imputed_self, na.rm = TRUE)
cleaned_pheno.data$edpgi_all_imputed_maternal <- (cleaned_pheno.data$edpgi_all_imputed_maternal
    - mean(cleaned_pheno.data$edpgi_all_imputed_self, na.rm = TRUE)) / sd(cleaned_pheno.data$edpgi_all_imputed_self, na.rm = TRUE)
cleaned_pheno.data$edpgi_exclude_imputed_paternal <- (cleaned_pheno.data$edpgi_exclude_imputed_paternal
    - mean(cleaned_pheno.data$edpgi_exclude_imputed_self, na.rm = TRUE)) / sd(cleaned_pheno.data$edpgi_exclude_imputed_self, na.rm = TRUE)
cleaned_pheno.data$edpgi_exclude_imputed_maternal <- (cleaned_pheno.data$edpgi_exclude_imputed_maternal
    - mean(cleaned_pheno.data$edpgi_exclude_imputed_self, na.rm = TRUE)) / sd(cleaned_pheno.data$edpgi_exclude_imputed_self, na.rm = TRUE)

# Scale Raw PGIs by the wider distribution. 
cleaned_pheno.data$asthma_pgi <- as.numeric(scale(cleaned_pheno.data$asthma_pgi))
cleaned_pheno.data$bipolar_pgi <- as.numeric(scale(cleaned_pheno.data$bipolar_pgi))
cleaned_pheno.data$bmi_pgi <- as.numeric(scale(cleaned_pheno.data$bmi_pgi))
cleaned_pheno.data$height_pgi <- as.numeric(scale(cleaned_pheno.data$height_pgi))
cleaned_pheno.data$schizophrenia_pgi <- as.numeric(scale(cleaned_pheno.data$schizophrenia_pgi))
cleaned_pheno.data$t2diabetes_pgi <- as.numeric(scale(cleaned_pheno.data$t2diabetes_pgi))
cleaned_pheno.data$adhd_pgi <- as.numeric(scale(cleaned_pheno.data$adhd_pgi))
cleaned_pheno.data$edpgi_all_raw_self <- as.numeric(scale(cleaned_pheno.data$edpgi_all_raw_self))
cleaned_pheno.data$edpgi_exclude_raw_self <- as.numeric(scale(cleaned_pheno.data$edpgi_exclude_raw_self))
cleaned_pheno.data$edpgi_all_imputed_self <- as.numeric(scale(cleaned_pheno.data$edpgi_all_imputed_self))
cleaned_pheno.data$edpgi_exclude_imputed_self <- as.numeric(scale(cleaned_pheno.data$edpgi_exclude_imputed_self))

#! Test: the genetic first-stage
lm(edpgi_all_imputed_self ~ 1 + I((edpgi_all_imputed_paternal + edpgi_all_imputed_maternal) / 2),
    data = cleaned_pheno.data) %>%
    summary() %>%
    print()
lm(edpgi_exclude_imputed_self ~ 1 + I((edpgi_exclude_imputed_paternal + edpgi_exclude_imputed_maternal) / 2),
    data = cleaned_pheno.data) %>%
    summary() %>%
    print()
lm(edpgi_all_imputed_self ~ 1 + I((edpgi_exclude_imputed_paternal + edpgi_exclude_imputed_maternal) / 2),
    data = cleaned_pheno.data) %>%
    summary() %>%
    print()
lm(edpgi_exclude_imputed_self ~ 1 + I((edpgi_all_imputed_paternal + edpgi_all_imputed_maternal) / 2),
    data = cleaned_pheno.data) %>%
    summary() %>%
    print()

# Show observation counts by both methods.
cleaned_pheno.data %>% filter(!is.na(edpgi_all_raw_self)) %>% NROW() %>% print()
cleaned_pheno.data %>% filter(!is.na(edpgi_all_imputed_self)) %>% NROW() %>% print()

# Validate the correlation between two different packages for computing Ed PGI.
print(cor(cleaned_pheno.data$edpgi_all_raw_self,
    cleaned_pheno.data$edpgi_all_imputed_self,
    use = "pairwise.complete.obs"))
print(cor(cleaned_pheno.data$edpgi_exclude_raw_self,
    cleaned_pheno.data$edpgi_exclude_imputed_self,
    use = "pairwise.complete.obs"))

#TODO: find a way to do both predict, feasible FEs, and na exclude.
# Fill in the hours work week variable, with imputed values.
#hours.reg <- lm(hours_workweek ~ 1 + sex_male
#    + edqual_highered                + edqual_professional           
#    + edqual_alevels                 + edqual_vocational             
#    + edqual_gcses                   + edqual_minimum                
#    + edqual_missing                 + edyears                       
#    + householdincome_less18k        + householdincome_18to31k
#    + householdincome_31to52k        + householdincome_52to100k
#    + householdincome_above100k      + householdincome_missing
#    | birthyear + jobcode_soc,
#    data = cleaned_pheno.data,
#    na.action = na.exclude)
#
#cleaned_pheno.data[is.na(cleaned_pheno.data$hours_workweek), ]
#table(round(predict(hours.reg, na.action = na.exclude,
#    newdata = cleaned_pheno.data[is.na(cleaned_pheno.data$hours_workweek), ])),
#    exclude = NULL)

## Fill in missing values row by row
#for (i in 1:nrow(cleaned_pheno.data)){
#    if (!is.na(cleaned_pheno.data[i, "hours_workweek"])){
#        next
#    }
#}


################################################################################
## Declare base data file.

# Rename and restrict to necessary columns.
final_pheno.data <- cleaned_pheno.data %>%
    mutate(
        # Ed PGI among entire sample
        edpgi_all_raw_self = edpgi_all_raw_self,
        edpgi_exclude_raw_self = edpgi_exclude_raw_self,
        # Ed PGI among sample with siblings -> main PGI in my analysis
        edpgi_all_imputed_self = edpgi_all_imputed_self,
        edpgi_all_imputed_paternal = edpgi_all_imputed_paternal,
        edpgi_all_imputed_maternal = edpgi_all_imputed_maternal,
        edpgi_exclude_imputed_self = edpgi_exclude_imputed_self,
        edpgi_exclude_imputed_paternal = edpgi_exclude_imputed_paternal,
        edpgi_exclude_imputed_maternal = edpgi_exclude_imputed_maternal,
        # Annual wage by hours worked, in thousands
        #TODO: make this vary by full/part-time
        soc_mean_annual = (soc_mean_hourly * hours_workweek * 40) / 1000) %>%
    # Mean parental Ed PGI, scaled by 
    rowwise() %>%
    mutate(edpgi_all_imputed_parental = mean(c(
            edpgi_all_imputed_paternal, edpgi_all_imputed_maternal), na.rm = TRUE),
        edpgi_exclude_imputed_parental = mean(c(
            edpgi_exclude_imputed_paternal, edpgi_exclude_imputed_maternal), na.rm = TRUE),
        ) %>%
    ungroup() %>%
    tibble() %>%
    # Mark the analysis sample, those with imputed Ed PGI + Ed + SOC data.
    mutate(analysis_sample = as.integer(!is.na(edpgi_all_imputed_parental)
        & !is.na(edyears) & !is.na(soc_mean_annual))) %>%
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
        # Ed PGI variables.
        edpgi_all_raw_self,
        edpgi_all_imputed_self,
        edpgi_all_imputed_parental,
        edpgi_all_imputed_paternal,
        edpgi_all_imputed_maternal,
        edpgi_exclude_raw_self,
        edpgi_exclude_imputed_self,
        edpgi_exclude_imputed_parental,
        edpgi_exclude_imputed_paternal,
        edpgi_exclude_imputed_maternal,
        # Genetic variables.
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
final_pheno.data %>% filter(!is.na(soc_mean_annual)) %>% NROW() %>% print()
final_pheno.data %>%
    filter(!is.na(edpgi_all_imputed_parental)
        & !is.na(edyears), !is.na(soc_mean_annual)) %>%
    NROW() %>%
    print()


################################################################################
## Calculate the random component of Ed PGI.

# Get the sample for whom we are getting mean expected PGI.
analysis.data <- final_pheno.data %>%
    filter(analysis_sample == 1)

# First the 23+me weights for the Ed PGI.
# (1) Calculate E[ Ed PGI | parents]
expected_edpgi_all.reg <- analysis.data %>%
    lm(edpgi_all_imputed_self ~ 1 + (
        poly(edpgi_all_imputed_paternal, 3) * poly(edpgi_all_imputed_maternal, 3)
        ) * sex_male * father_present * mother_present,
        na.action = na.exclude,
        data = .)
# (2) Calculate random component = Ed PGI - \hat E[ Ed PGI | parents]
analysis.data$edpgi_all_imputed_random <- (analysis.data$edpgi_all_imputed_self
    - predict(expected_edpgi_all.reg, analysis.data))

# Second the UKB weights for the Ed PGI.
# (1) Calculate E[ Ed PGI | parents]
expected_edpgi_exclude.reg <- analysis.data %>%
    lm(edpgi_exclude_imputed_self ~ 1 + (
        poly(edpgi_exclude_imputed_paternal, 3) * poly(edpgi_exclude_imputed_maternal, 3)
        ) * sex_male * father_present * mother_present,
        na.action = na.exclude,
        data = .)
# (2) Calculate random component = Ed PGI - \hat E[ Ed PGI | parents]
analysis.data$edpgi_exclude_imputed_random <- (analysis.data$edpgi_exclude_imputed_self
    - predict(expected_edpgi_exclude.reg, analysis.data))


#! Testing the first-stage ORIV
analysis.data %>%
    lm(edpgi_exclude_imputed_self ~ 1 + edpgi_all_imputed_random,
        na.action = na.exclude,
        data = .) %>%
    summary() %>%
    print()

# Show summary statistics of the random component.
print(mean(analysis.data$edpgi_all_imputed_random))
print(sd(analysis.data$edpgi_all_imputed_random))
hist(analysis.data$edpgi_all_imputed_random)

# Put Ed PGi random component back onto the phenotype data file.
final_pheno.data <- analysis.data %>%
    select(eid, edpgi_all_imputed_random, edpgi_exclude_imputed_random) %>%
    right_join(final_pheno.data, by = "eid")


################################################################################
## Connect UKB data with collected distance to nearest uni data.

# Read the higher ed locations file.
higher_loc.data <- read_csv(
    file.path(input.folder, "..", "..", "uk-highered", "highered-compiled.csv"))

# Get the subsample with birth locations.
uni.age <- 16
birth.data <- final_pheno.data %>%
    filter(!is.na(birth_n_coord), !is.na(birth_e_coord),
        !is.na(birthyear), analysis_sample == 1) %>%
    mutate(year_ageuni = birthyear + uni.age) %>%
    select(eid, birth_n_coord, birth_e_coord, year_ageuni)
# Empty columns for uni data
# Uni if open at age when deciding uni.
birth.data$open_closest_uni_name <- NA
birth.data$open_closest_uni_county <- NA
birth.data$open_year_founded <- NA
birth.data$open_uni_n_coord <- NA
birth.data$open_uni_e_coord <- NA
# all unis
birth.data$all_closest_uni_name <- NA
birth.data$all_closest_uni_county <- NA
birth.data$all_year_founded <- NA
birth.data$all_uni_n_coord <- NA
birth.data$all_uni_e_coord <- NA

# Loop across each individual.
total.rows <- nrow(birth.data)
for (i in 1:total.rows) {
    if (((100 * i / round(total.rows, -2)) %% 1) == 0) {
        print(paste0(i, " out of ", total.rows, ", ", 100 * (i / total.rows), "% done."))
    }
    #print(paste0(i, " out of ", total.rows, ", ", 100 * (i / total.rows), "% done."))
    individual.data <- birth.data[i, ]
    # Get universities that existed when this person was 18
    open_unis.data <- higher_loc.data[
        higher_loc.data$year_founded <= individual.data$year_ageuni, ]
    # Find which open uni is the closest to individual i 
    open_uni_dist.list <- sqrt(
        (open_unis.data$uni_n_coord - individual.data$birth_n_coord)^2 + (
            open_unis.data$uni_e_coord - individual.data$birth_e_coord)^2)
    closest_open_uni.index <- which.min(open_uni_dist.list)
    # Save this info to the birth data file.
    closest_open_uni.data <- open_unis.data[closest_open_uni.index, ]
    birth.data[i, ]$open_closest_uni_name <- closest_open_uni.data$uni_name
    birth.data[i, ]$open_closest_uni_county <- closest_open_uni.data$uni_county_name
    birth.data[i, ]$open_year_founded <- closest_open_uni.data$year_founded
    birth.data[i, ]$open_uni_n_coord <- closest_open_uni.data$uni_n_coord
    birth.data[i, ]$open_uni_e_coord <- closest_open_uni.data$uni_e_coord
    # Do the same for all universities (even if not open yet).
    uni_dist.list <- sqrt(
        (higher_loc.data$uni_n_coord - individual.data$birth_n_coord)^2 + (
            higher_loc.data$uni_e_coord - individual.data$birth_e_coord)^2)
    all_uni.index <- which.min(uni_dist.list)
    # Save this info to the birth data file.
    closest_all_uni.data <- higher_loc.data[all_uni.index, ]
    birth.data[i, ]$all_closest_uni_name <- closest_all_uni.data$uni_name
    birth.data[i, ]$all_closest_uni_county <- closest_all_uni.data$uni_county_name
    birth.data[i, ]$all_year_founded <- closest_all_uni.data$year_founded
    birth.data[i, ]$all_uni_n_coord <- closest_all_uni.data$uni_n_coord
    birth.data[i, ]$all_uni_e_coord <- closest_all_uni.data$uni_e_coord
}


################################################################################
## Add counties to birth places.

# Read UK district (i.e., county) borders in 1971, from UK data service.
# https://borders.ukdataservice.ac.uk/bds.html
eng_shape.data <- file.path(input.folder, "..", "..", "uk-highered",
        "locations", "England_ct_1971", "england_ct_1971.shp") %>% 
    st_read() %>%
    mutate(country = "England")
scot_shape.data <- file.path(input.folder, "..", "..", "uk-highered",
        "locations", "Scotland_ct_1971", "scotland_ct_1971.shp") %>% 
    st_read()  %>%
    mutate(country = "Scotland")
wales_shape.data <- file.path(input.folder, "..", "..", "uk-highered",
        "locations", "Wales_ct_1971", "wales_ct_1971.shp") %>% 
    st_read() %>%
    mutate(country = "Wales")
# Combine, to get entire UK data.
uk_shape.data <- rbind(eng_shape.data, scot_shape.data, wales_shape.data)

# Ensure both datasets use the same coordinate reference system
print(paste("UKB birth data CRS:", st_crs(birth.data)$input))
print(paste("County data CRS:", st_crs(uk_shape.data)$input))

# Put the birth county onto the birth data.
birth.data <- birth.data %>%
    st_as_sf(coords = c("birth_e_coord", "birth_n_coord"),
        crs = st_crs(uk_shape.data)) %>%
    st_join(uk_shape.data, join = st_within)

# Restrict to relevant columns
birth.data <- birth.data %>%
    rename(
        #birth_county_label = label,
        birth_county_name = name) %>%
    tibble() %>%
    select(-country, -label, -geometry, -year_ageuni)

# Put the ed location data back on to the dataframe.
final_pheno.data <- final_pheno.data %>%
    left_join(birth.data, by = "eid")

# How many have newly opened unis?
final_pheno.data %>%
    filter(analysis_sample == 1) %>%
    NROW() %>%
    print()
final_pheno.data %>% 
    filter(analysis_sample == 1) %>%
    filter(open_closest_uni_name != all_closest_uni_name) %>%
    NROW() %>%
    print()


################################################################################
## Save the resulting data file.

# Save the file.
final_pheno.data %>%
    write_csv(file.path(output.folder, "ukb-cleaned-pheno.csv"))
