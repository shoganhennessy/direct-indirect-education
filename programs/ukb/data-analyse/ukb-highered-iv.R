#!/usr/bin/R
## Senan Hogan-Hennessy, 29 May 2025
## UKB data -> First-stage analysis of uni openings/distance -> education.
set.seed(47)
print(Sys.time())

# Functions for data manipulation and visualisation
library(tidyverse)
# Functions for tables into TeX
library(xtable)
# The standard, linear, IV estimator package.
library(fixest)
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

#! Final test -> get population counts for the county born in., in year choosing uni)
analysis.data$year_ageuni <- (analysis.data$birthyear + 17) #1955--1985
hist(analysis.data$year_ageuni)
table(analysis.data$birth_county_name, exclude = NULL)
population.data <- read_csv(
    file.path(data.folder, "..", "..", "uk-highered", "locations",
        "Counties_with_1951_population_merged.csv")) %>%
    transmute(
        birth_county_name = county,
        birth_county_pop1951 = pop1951,
        birth_county_pop1961 = pop1961,
        birth_county_pop1971 = pop1971,
        birth_county_pop1981 = pop1981,
        birth_county_pop1985 = pop1981) %>% 
    pivot_longer(
        cols = starts_with("birth_county_pop"),
        names_to = "year",
        names_prefix = "birth_county_pop",
        values_to = "birth_county_pop",
        values_drop_na = FALSE) %>%
    mutate(year = as.integer(year))
# Then expand across years, can change as needed in future.
population.data <- population.data %>%
    expand(birth_county_name, year = 1951:1985) %>%
    left_join(population.data, by = c("birth_county_name", "year")) %>%
    fill(birth_county_pop, .direction = "down")
# Add on to the analysis data
analysis.data <- analysis.data %>%
    left_join(population.data, join_by(
        birth_county_name == birth_county_name,
        year_ageuni == year))

# Counting the number of unis in your birth county.
higher_loc.data <- read_csv(
    file.path(data.folder, "..", "..", "uk-highered", "highered-compiled.csv"))

# Count new unis over the right time period.
yearly_intake.data <- higher_loc.data %>%
    filter(year_founded <= 1950) %>% 
    group_by(uni_county_name, year_founded) %>% 
    summarise(new_unis = n(), .groups = "drop")
all_years <- 1000:1970
county_year.data <- expand_grid(
    uni_county_name = sort(unique(yearly_intake.data$uni_county_name)),
    year_founded = all_years)
higher_count.data <- county_year.data %>% 
    left_join(yearly_intake.data, by = c("uni_county_name", "year_founded")) %>% 
    mutate(new_unis = replace_na(new_unis, 0L)) %>%
    arrange(uni_county_name, year_founded) %>% 
    group_by(uni_county_name) %>% 
    mutate(open_uni_count = cumsum(new_unis)) %>% 
    ungroup() %>%
    select(uni_county_name, year_founded, open_uni_count)
# Join to the analysis.data
analysis.data <- analysis.data %>%
    left_join(higher_count.data, by = join_by(
        birth_county_name == uni_county_name, birthyear == year_founded)) %>%
    mutate(open_uni_count = replace_na(open_uni_count, 0L))

# Clean the open uni per person data.
analysis.data <- analysis.data %>%
    mutate(open_uni_perperson = open_uni_count / (birth_county_pop / 10^3))
# Upper outlier.
upper.limit <- quantile(analysis.data$open_uni_perperson, 0.99, na.rm = TRUE)
analysis.data <- analysis.data %>%
    mutate(open_uni_perperson = ifelse(open_uni_perperson > upper.limit,
        NA, open_uni_perperson))
hist(analysis.data$open_uni_perperson)
mean(analysis.data$open_uni_perperson, na.rm = TRUE)


# Calculate distance, and ed rate in your birthcounty.
bin.size <- 5
analysis.data <- analysis.data %>%
    mutate(
        open_uni_birth_county =  as.integer(
            open_closest_uni_county != birth_county_name),
        open_uni_birth_dist = sqrt(
            (open_uni_n_coord - birth_n_coord)^2 + (
                open_uni_e_coord - birth_e_coord)^2) / 1000,
        all_uni_birth_county = as.integer(
            all_closest_uni_county != birth_county_name),
        all_uni_birth_dist = sqrt(
            (all_uni_n_coord - birth_n_coord)^2 + (
                all_uni_e_coord - birth_e_coord)^2) / 1000) %>%
    mutate(birthyear_bin = bin.size * round(birthyear / bin.size)) %>%
    group_by(birth_county_name, birthyear_bin) %>%
    mutate(highered_rate_county = mean(edqual_highered * 100)) %>%
    ungroup()

# Test the first stage.
analysis.data %>%
    feols(edyears ~ 1 
        + highered_rate_county
        + open_uni_perperson
        + open_uni_birth_dist
        + birth_county_pop
        | birthyear,
        data = .) %>%
    summary() %>%
    print()
# Implied IV
analysis.data %>%
    feols(log(soc_mean_annual) ~ 1
        + birth_county_pop
        | birthyear
        | edyears ~ highered_rate_county + open_uni_perperson + open_uni_birth_dist + birth_county_pop,
        data = .) %>%
    summary() %>%
    print()
# Just for higher ed.
analysis.data %>%
    filter(14 <= edyears, edyears <= 18) %>%
    feols(edqual_highered ~ 1 + birth_county_pop
        + highered_rate_county + open_uni_perperson + open_uni_birth_dist
        | birthyear + birth_county_name,
        data = .) %>%
    summary() %>%
    print()
# Implied IV
analysis.data %>%
    filter(14 <= edyears, edyears <= 18) %>%
    feols(log(soc_mean_hourly) ~ 1 + birth_county_pop
        | birthyear + birth_county_name
        | edpgi_all_imputed_self + edyears ~
            edpgi_all_imputed_random + highered_rate_county,
        data = .) %>%
    summary() %>%
    print()


analysis.data %>%
    filter(edyears <= 14) %>%
    mutate(highschool = as.integer(14 <= edyears)) %>%
    feols(log(soc_mean_hourly) ~ 1 + birth_county_pop
        | birthyear
        | edyears ~ highered_rate_county + open_uni_perperson + open_uni_birth_dist,
        data = .) %>%
    summary() %>%
    print()



################################################################################
## Sibling difference returns to education.

# Get a control for youngest.
analysis.data <- analysis.data %>%
    group_by(famid) %>%
    mutate(youngest_sibling =
        as.integer(birthyear + birthmonth / 12
            == min(birthyear + birthmonth / 12)),
        edmin = min(edyears),
        edmax = max(edyears),
        eddiff = edmax - edmin) %>%
    ungroup()



# Effect of PGI on higher ed (first-stage)
analysis.data %>%
    #filter(edmin != edmax) %>%
    feols(edyears ~ 1
        + adhd_pgi + asthma_pgi + bipolar_pgi + bmi_pgi + height_pgi
        + schizophrenia_pgi + t2diabetes_pgi + sex_male
        + factor(visityear) + poly(recruitedage, 3) + factor(sibling_count) + urban
        | edpgi_all_imputed_self ~ edpgi_exclude_imputed_random,
        data = .) %>%
    summary() %>%
    print()

source("oriv.R")
oriv.iv <- GORIV(edyears ~ 1
        + adhd_pgi + asthma_pgi + bipolar_pgi + bmi_pgi + height_pgi
        + schizophrenia_pgi + t2diabetes_pgi + sex_male
        + factor(visityear) + poly(recruitedage, 3) + factor(sibling_count) + urban,
    "edpgi_all_imputed_self", "edpgi_exclude_imputed_self",
    IV.list = c(
        "edpgi_all_imputed_random", "edpgi_exclude_imputed_random"),
    control.list = c(
        "edpgi_all_imputed_parental", "edpgi_exclude_imputed_parental"),
    IID = "eid", FID = "famid",
    data = analysis.data)
print(summary(oriv.iv))

# Effect of higher ed on income (second-stage)
analysis.data %>%
    filter(edmin != edmax) %>%
    feols(log(soc_mean_hourly) ~ 1
        + edyears
        + edpgi_all_imputed_parental
        | famid + birthyear + sex_male^youngest_sibling
        | edpgi_all_imputed_self ~ edpgi_exclude_imputed_self,
        data = .) %>%
    summary() %>%
    print()

# How many families have education differences?
sibdiff.data <- analysis.data %>%
    group_by(famid) %>%
    summarise(
        edmin = min(edyears),
        edmax = max(edyears)) %>%
    ungroup()
sibdiff.data
sibdiff.data %>% filter(edmin != edmax)
sibdiff.data %>% filter(edmin != edmax, edmax == 18)


################################################################################
## Testing out the first-stage, with binary measure of uni in local labour market.

# 1. Binary for uni in that county when age 18.
mean(analysis.data$all_uni_birth_county, na.rm = TRUE)
mean(analysis.data$open_uni_birth_county, na.rm = TRUE)
analysis.data %>%
    feols(edyears ~ 1 + open_uni_birth_county * all_uni_birth_county
        + birth_county_pop
        | birthyear + birth_county_name,
        data = .) %>%
    summary() %>%
    print()
# IV
mean(analysis.data[analysis.data$edqual_highered == 0, ]$soc_mean_hourly)
analysis.data %>%
    feols(log(soc_mean_hourly) ~ 1
        | birthyear
        | edyears ~ open_uni_birth_county,
        data = .) %>%
    summary() %>%
    print()


# 2. Distance to the nearest open university.
distance.reg <- analysis.data %>%
    feols(edqual_highered ~ 1 +
        log(open_uni_birth_dist)
        + birth_county_pop
        | birthyear + birth_county_name,
        data = .)
distance.reg %>%
    summary() %>%
    print()

# 3. Years until uni opening locally (quasi-event study).
analysis.data <- analysis.data %>%
    mutate(uni_opening = (birthyear + 17) - all_year_founded)

# Quasi-DiD for university opening.
# DiD data, around the university opening in the local district.
max.years <- 11
did.data <- analysis.data %>%
    filter(all_uni_birth_county == 0 | (
        all_uni_birth_county == 1 & abs(uni_opening) < max.years)) %>%
    mutate(uni_opening = ifelse(all_uni_birth_county == 0, 0, uni_opening))

did.data %>% filter(-max.years <= uni_opening & uni_opening <= 0) %>% nrow()
did.data %>% filter(0 < uni_opening & uni_opening <= max.years) %>% nrow()
# Event- study visualisation of it.
did.data %>%
    lm(edqual_highered ~ 1 +
        factor(uni_opening, levels = c(0, -max.years:-1, 1:max.years)
            ),
        data = .) %>%
    summary() %>%
    print()
# DiD specification, too.
did.data %>%
    feols(edyears ~ 1 + I(uni_opening > 0) * all_uni_birth_county 
        | birthyear,
        data = .) %>%
    summary() %>%
    print()
# Implied returns to higher education.
library(fixest)
did.data %>%
    feols(log(soc_mean_hourly) ~ 1 +
        I(uni_opening > 0) + all_uni_birth_county 
        | birthyear
        | edqual_highered ~ I(uni_opening > 0) : all_uni_birth_county ,
        data = .) %>%
    summary() %>%
    print()

# 4. Continuous measure of distance from birth to open uni.
analysis.data %>%
    feols(edqual_highered ~ 1 +
        + log(open_uni_birth_dist)
        | birthyear,
        data = .) %>%
    summary() %>%
    print()
# Implied IV
library(fixest)
analysis.data %>%
    feols(log(soc_mean_hourly) ~ 1
        | birthyear + birth_county_name
        | edqual_highered ~ open_uni_birth_county,
        data = .) %>%
    summary() %>%
    print()

#! Testing it out: counting the number of unis in your county.
higher_loc.data <- read_csv(
    file.path(data.folder, "..", "..", "uk-highered", "highered-compiled.csv"))


# Count new unis over the right time period.
yearly_intake.data <- higher_loc.data %>%
    filter(year_founded <= 1950) %>% 
    group_by(uni_county_name, year_founded) %>% 
    summarise(new_unis = n(), .groups = "drop")
all_years <- 1000:1970
county_year.data <- expand_grid(
    uni_county_name = sort(unique(yearly_intake.data$uni_county_name)),
    year_founded = all_years)
higher_count.data <- county_year.data %>% 
    left_join(yearly_intake.data, by = c("uni_county_name", "year_founded")) %>% 
    mutate(new_unis = replace_na(new_unis, 0L)) %>%
    arrange(uni_county_name, year_founded) %>% 
    group_by(uni_county_name) %>% 
    mutate(open_uni_count = cumsum(new_unis)) %>% 
    ungroup() %>%
    select(uni_county_name, year_founded, open_uni_count)
# Join to the analysis.data
final_pheno.data <- final_pheno.data %>%
    left_join(higher_count.data, by = join_by(
        birth_county_name == uni_county_name, birthyear == year_founded)) %>%
    mutate(open_uni_count = replace_na(open_uni_count, 0L))

# Clean the open uni per person data.
final_pheno.data <- final_pheno.data %>%
    mutate(open_uni_perperson = (open_uni_count / (birth_county_pop / 100))) %>%
    mutate(open_uni_perperson = ifelse((open_uni_count / birth_county_pop) > 0.1, NA,
        open_uni_perperson))
    
# Test the first stage.
final_pheno.data %>%
    filter(analysis_sample == 1) %>%
    feols(edyears ~ 1
        + open_uni_perperson
        | birth_county_name + birthyear,
        data = .) %>%
    summary() %>%
    print()
# Implied IV
library(fixest)
final_pheno.data %>%
    filter(analysis_sample == 1) %>%
    feols(log(soc_mean_hourly) ~ 1
        | birth_county_name
        | edyears ~ I(open_uni_count),
        data = .) %>%
    summary() %>%
    print()
final_pheno.data %>%
    filter(analysis_sample == 1, 14 <= edyears, edyears <= 18) %>%
    feols(log(soc_mean_hourly) ~ 1
        | birthyear + birth_county_name
        | edqual_highered ~ I(open_uni_count / birth_county_pop),
        data = .) %>%
    summary() %>%
    print()

#TODO Findings:
#TODO People born in a county with a uni have less education.  Inner city worse off than rural.
#TODO Last ditch effort, get data on county/district population, and use the division as the IV.
#TODO This is the same as Currie (QJE 2004).
