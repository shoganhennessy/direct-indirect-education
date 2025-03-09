library(data.table)

# READ data (tsv or csv format)
UKB <- 
UKB <- data.table(UKB)

############################################################################################3
############################################################################################3
### required info

UKB[, male := f.31.0.0] 
UKB[, age := as.numeric(f.21003.0.0)]
UKB[, year := substr(f.53.0.0, 1, 4)]  # year of assessment
UKB[, yob := as.numeric(f.34.0.0)]

############################################################################################3
############################################################################################3
### Curent job
UKB[, SOC := f.20277.0.0] # the field may be different, depending on the version of your data
UKB[, SOC := substr(f.20277.0.0, 1, 4)] 

        # if older version of UKB data (132, 20024 present)
        # UKB[, SOC := f.132.0.0] # the field may be different, depending on the version of your data
        # UKB[, SOC := substr(f.132.0.0, 1, 4)] 
        # UKB[, SOC := ifelse(SOC==0, f.20024.0.0, SOC)]


############################################################################################3
############################################################################################3
### Get SOC from past job history
job_temp <- as.matrix(UKB[, .SD, .SDcols=grep("f.22617.0.", names(UKB), value=T)]) 
year_s <- as.matrix(UKB[, .SD, .SDcols=grep("f.22602.0.", names(UKB), value=T)]) # year job started
year_e <- as.matrix(UKB[, .SD, .SDcols=grep("f.22603.0.", names(UKB), value=T)]) # year job ended

year_s[year_s<0] <- NA
year_e[year_e<0] <- NA

# year_e missing if currently held job. Give 2010 if job started before 2010
year_e <- matrix(ifelse(!is.na(job_temp) & is.na(year_e) & year_s <=2010, 2010, year_e), ncol=40)

# year_e cap to 2010
year_e <- matrix(ifelse(year_e > 2010, 2010, year_e), ncol=40)

age_s <- year_s - UKB$yob
age_e <- year_e - UKB$yob

# leave only job held between age 40 and 64. 
job_temp <- matrix(ifelse((age_s < 40 & age_e %between% c(40,64)) | (age_s %in% 40:64 & age_e <=64), job_temp, NA), ncol=40)

# leave only job held between year 2002 and 2010
job_temp <- matrix(ifelse( (year_s<2002 & year_e>=2002) | (year_s %between% c(2002,2010)), job_temp, NA), ncol=40)

# get column index to extract (latest job info available)
ind <- apply(job_temp, 1, function(x) max(which(!is.na(x)), na.rm=T) )
ind <- ifelse(ind==-Inf, NA, ind)

UKB$past_SOC <-  sapply(1:nrow(job_temp), function(x) ifelse(!is.na(ind[x]), job_temp[x,ind[x]], NA))
UKB$past_SOC_age <- sapply(1:nrow(age_e), function(x) ifelse(!is.na(ind[x]), age_e[x,ind[x]], NA))
UKB$past_SOC_year <- sapply(1:nrow(year_e), function(x) ifelse(!is.na(ind[x]), year_e[x,ind[x]], NA))



############################################################################################3
############################################################################################3
### Use past SOC for those too old (>64) as of 2002-2010, or current SOC missing
UKB[, SOC := ifelse(age > 64, NA, SOC)]

UKB[, use_past_SOC := ifelse(is.na(SOC) & !is.na(past_SOC), 1, 0) ]
UKB[, SOC := ifelse(use_past_SOC==1, past_SOC, SOC)]
# change age, year to match past SOC
UKB[, age := ifelse(use_past_SOC==1, past_SOC_age, age)]
UKB[, year := ifelse(use_past_SOC==1, past_SOC_year, year)]



############################################################################################3
############################################################################################3
### Get occupational wages

source("./function_script.R")

UKB <- GET_INCOME(UKB)