#!/usr/bin/R
## Senan Hogan-Hennessy, 5 AUgust 2024.
## Script to ingest raw HRS + HRS/RAND data, output usable format.
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

# Load the HRS data panel (previously compiled).
hrs_panel.data <- output.folder %>%
    file.path("hrs-panel-data.csv") %>%
    read_csv()

# Load the HRS gene score data, for those with Euro ancestory
gene_euro.data <- data.folder %>%
    file.path("hrs-public", "PGENSCORE4r3", "pgenscore4e_r.dta") %>%
    read_dta()

# Load the early childhood supplement.
childhood.data <- data.folder %>%
    file.path("hrs-public",
        "ChildhoodHealthAndFamily", "AAGGCHLDFH2016A_R.dta") %>%
    read_dta()

# CPI-U from FREDS. Yearly average, seasonally adjusted, base 1982-1984=100
cpiu.data <- data.folder %>%
    file.path("inflation-indices", "freds-cpiu.csv") %>%
    read_csv()


################################################################################
## Clean early childhood factors data.

# Collect the relevant variables from childhood SES supplement.
clean_childhood.data <- childhood.data %>%
    # Note: VAR %in% c(list) syntax means NAs in VAR go to 0 (not NA)
    transmute(
        # Household + person identifiers. 
        rahhid = HHID,
        rahhidpn = paste0(HHID, PN),
        # Rating of family income
        family_poor = ifelse(FAMFIN %in% c(5, 6), 1, 0),
        # Move to due to fin difficulties
        family_move = ifelse(MOVFIN %in% c(1), 1, 0),
        # Family receive financial help (from relatives) during childhood.
        family_finhelp = ifelse(FMFINH %in% c(1), 1, 0),
        # Father missing
        father_missing = ifelse(FAUNEM %in% c(7), 1, 0),
        # Father unemployed during childhood.
        father_unemp = ifelse(FAUNEM %in% c(1), 1, 0),
        # Father occupation code
        father_manualjob = ifelse(FJOB %in% c(5), 1, 0),
        # Father + mother ed years
        father_edyears_supp = ifelse(FAEDUC %in% c(0:17), FAEDUC, NA),
        mother_edyears_supp = ifelse(MOEDUC %in% c(0:17), MOEDUC, NA),
        # Child health rating
        child_badhealth = ifelse(RTHLTHCH %in% c(5, 6), 1, 0),
        # Parents smoke
        parents_smoke = ifelse(PARSMOKE %in% c(1), 1, 0),
        # Childhood head injury
        child_headinjury = ifelse(CHHEADINJ %in% c(1), 1, 0),
        # CHildhood disability
        child_disability = ifelse(CHDISABL %in% c(1), 1, 0))


################################################################################
## Clean genetic data.

# Standardise the HRS gene data, with values for Euro ancestory
clean_gene.data <- gene_euro.data %>%
    # Select the relevant columns
    transmute(
        # Standardise the unique identifying number.
        rahhidpn = paste0(hhid, pn),
        # Gene score Lee+(2018) educ years.
        genescore_educ_euro = E4_EA3_W23_SSGAC18,
        # Gene score Davies+(2018) general cognition.
        genescore_gencog_euro = E4_GCOG2_CHARGE18,
        # Gene score GIANT(2018) BMI.
        genescore_bmi_euro = E4_BMI2_GIANT18,
        # Gene score GIANT(2015) waist-to-hips ratio.
        genescore_waisthips_euro = E4_WHR_GIANT15,
        # Gene score GIANT(2018) Height.
        genescore_height_euro = E4_HEIGHT2_GIANT18,
        # Gene score PGC(2014) Schizophrenia.
        genescore_schizo_euro = E4_SCZ_PGC14,
        # Gene score PGC(2018) Alcohol dependence.
        genescore_alcohol_euro = E4_ALC_PGC18,
        # Gene score IGAP(2013) Alzheimers.
        genescore_alzh_euro = E4_AD_IGAP13,
        # Gene score SSGAC(2016) Neuroticism.
        genescore_neurotic_euro = E4_NEUROT_SSGAC16,
        # Gene score SSGAC(2016) Wellbeing.
        genescore_wellbeing_euro = E4_WELLB_SSGAC16,
        # Gene score CARDIoGRAM(2011) Coronary Artery Disease.
        genescore_heartcvd_euro = E4_CD_CARDIOGRAM11,
        # Gene score CARDIoGRAM(2011) Myocardial Infarction.
        genescore_heartmi_euro = E4_MI_CARDIOGRAM15,
        # Gene score PGC(2011) Bipolar Disorder.
        genescore_bipolar_euro = E4_BIP_PGC11,
        # Gene score PGC(2017) ADHD.
        genescore_adhd_euro = E4_ADHD_PGC17,
        # Gene score PGC(2013) Mental health cross disorder.
        genescore_crossdisorder_euro = E4_XDISORDER_PGC13,
        # Gene score GPC(2016) Extraversion Polygenic Score.
        genescore_extravert_euro = E4_EXTRAVER_GPC16,
        # Gene score CHARGE(2015) Longevity.
        genescore_longevity_euro = E4_LONG_CHARGE15,
        # Gene score BROAD(2017) Antisocial behaviour.
        genescore_antisocial_euro = E4_AB_BROAD17,
        # Gene score IOCDF(2017) Obsessive Compulsive Disorder.
        genescore_ocd_euro = E4_OCD_IOCDF17,
        # Gene score PGC(2018) Major Depressive Disorder 2.
        genescore_mdepressived_euro = E4_MDD2_PGC18,
        # Gene score Duncan+(2018) PTSD.
        genescore_ptsd_euro = E4_PTSDEA_PGC18,
        # Gene score Otowa+(2016) Anxiety factor score.
        genescore_anxietyfactor_euro = E4_ANXFS_ANGST16,
        # Gene score Otowa+(2016) Anxiety case control score.
        genescore_anxietycontrol_euro = E4_ANXCC_ANGST16,
        # Gene score GSCAN(2019) Age at smoking.
        genescore_agesmoking_euro = E4_AI_GSCAN19)


################################################################################
## Add on inflation index to the panel data.

# Clean CPI data, and put in terms of base year 2022.
base.year <- 2022
base.index <- cpiu.data %>%
    filter(base.year ==
        (cpiu.data$DATE %>% substr(start = 0, stop = 4) %>% as.integer())) %>%
    pull(CPIAUCSL)
cpi.data <- cpiu.data %>%
    transmute(
        survey_year = DATE %>% substr(start = 0, stop = 4) %>% as.integer(),
        cpi = CPIAUCSL / base.index)

# Connect the CPI to the panel data.
hrs_panel.data <- hrs_panel.data %>%
    left_join(cpi.data, by = "survey_year")


################################################################################
## Save the resulting data file.

# Save the early childhood data file.
clean_childhood.data %>%
    write_csv(file.path(output.folder, "hrs-childhood-data.csv"))

# Save the gene score data file.
clean_gene.data %>%
    write_csv(file.path(output.folder, "hrs-gene-data.csv"))

# Save the HRS panel data file.
hrs_panel.data %>%
    write_csv(file.path(output.folder, "hrs-panel-cleaned.csv"))

# Save the labels for the data file.
sink(file = file.path(data.folder, "hrs-panel-labels.txt"))
print("Here are the column labels for the HRS data:")
hrs_panel.data %>%
    select(starts_with("ra"),
        starts_with("r13"), starts_with("h13")) %>%
    sapply(attr, "label") %>%
    print()

print("Here are the labels for encoded values for the HRS data:")
hrs_panel.data %>%
    select(starts_with("ra"),
        starts_with("r13"), starts_with("h13")) %>%
    sapply(attr, "labels") %>%
    print()
sink(file = NULL)
