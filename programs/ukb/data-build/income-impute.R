#!/usr/bin/R
## Senan Hogan-Hennessy, 18 Feb 2025
## Convert the Kweon+ (2024) UKB income imputation file to a usable format.
# Source: 1. https://doi.org/10.62891/aac85602
#         2. https://doi.org/10.1038/s41562-024-02080-7


################################################################################
## Senan note: Code that follows is provided by Kweon+ (2024) rep package.
## Senan note: Only changes are making the indentation consistent.


# PLEASE ensure you have "data_input.RData" 

## Function to impute income from occupation code
# INPUT: data set (data.frame or data.table object) for which to impute income

# IMPORTANT : The data you will provide to the function must contain the following variables with the same name and class specified 
    # "SOC" (numeric or factor): 4-digit standardised occupation code (2000)
    # "male" (numeric) : dummy variable for being male
    # "age" (numeric or factor) : age, or (year of birth - year of observation); restricted to 35~64
    # "year" (numeric, character, or factor) : year of observation, format must be 'yyyy', only works for 2002 ~ 2016

# OUTPUT: data.table object with log occupational income (data.table object; the name of the log hourly occupational income is log_y_hourly)

# note: you need to have data.table package installed

GET_INCOME <- function(D){
    library(data.table)
    D <- data.table(D)
    D[, year:=factor(year)]
    y <- as.numeric(levels(D$year))
    for (i in y){
        ONS_temp <- as.data.frame(ONSlist[i-2001])  
        soc <- as.matrix(D[year==i, "SOC"])

        ###get mean of each occupation (by gender)
        rown <- match(soc, ONS_temp$code) 
        D[year==i,
            soc_mean_hourly:=ifelse(male == 1, ONS_temp[rown, "mean_male"],
                ONS_temp[rown, "mean_female"])]
        
        #if 4-digit info missing
        rown <- match(substr(soc, 1, 3), ONS_temp$code) 
        D[year==i,
            soc_mean_hourly:=ifelse(!is.na(soc) & is.na(soc_mean_hourly),
                ifelse(male == 1, ONS_temp[rown, "mean_male"],
                    ONS_temp[rown, "mean_female"]), soc_mean_hourly) ]

        # if 3-digit info also missing
        rown <- match(substr(soc, 1, 2), ONS_temp$code) 
        D[year==i,
            soc_mean_hourly:=ifelse(!is.na(soc) & is.na(soc_mean_hourly),
                ifelse(male == 1, ONS_temp[rown, "mean_male"],
                    ONS_temp[rown, "mean_female"]),
                        soc_mean_hourly) ]
        ###get median of each occupation (by gender)
        rown <- match(soc, ONS_temp$code) 
        D[year==i,
            soc_median_hourly:=ifelse(male == 1, ONS_temp[rown, "median_male"],
                ONS_temp[rown, "median_female"])]
        #if 4-digit info missing
        rown <- match(substr(soc, 1, 3), ONS_temp$code) 
        D[year==i,
            soc_median_hourly:=ifelse(!is.na(soc) & is.na(soc_median_hourly),
                ifelse(male == 1, ONS_temp[rown, "median_male"],
                    ONS_temp[rown, "median_female"]),
                        soc_median_hourly) ]
        # if 3-digit info also missing
        rown <- match(substr(soc, 1, 2), ONS_temp$code) 
        D[year==i,
            soc_median_hourly:=ifelse(!is.na(soc) & is.na(soc_median_hourly),
                ifelse(male == 1, ONS_temp[rown, "median_male"],
                    ONS_temp[rown, "median_female"]),
                        soc_median_hourly) ]
        ### get mean for all
        rown <- match(soc, ONS_temp$code)
        D[year==i, soc_mean_hourly_all:= ONS_temp[rown, "mean_all"]]

        #if 4-digit info missing
        rown <- match(substr(soc, 1, 3), ONS_temp$code) 
        D[year==i,
            soc_mean_hourly_all:= ifelse(!is.na(soc) &
                is.na(soc_mean_hourly_all), ONS_temp[rown, "mean_all"],
                    soc_mean_hourly_all) ]
        # if 3-digit info also missing
        rown <- match(substr(soc, 1, 2), ONS_temp$code)
        D[year==i,
            soc_mean_hourly_all:= ifelse(!is.na(soc) &
                is.na(soc_mean_hourly_all), ONS_temp[rown, "mean_all"],
                    soc_mean_hourly_all) ]
        ### get median for all
        rown <- match(soc, ONS_temp$code) 
        D[year==i, soc_median_hourly_all:= ONS_temp[rown, "median_all"]]

        #if 4-digit info missing
        rown <- match(substr(soc, 1, 3), ONS_temp$code) 
        D[year==i,
            soc_median_hourly_all:= ifelse(!is.na(soc) &
                is.na(soc_median_hourly_all), ONS_temp[rown, "median_all"],
                    soc_median_hourly_all) ]
        
        # if 3-digit info also missing
        rown <- match(substr(soc, 1, 2), ONS_temp$code) 
        D[year==i,
            soc_median_hourly_all:= ifelse(!is.na(soc) &
                is.na(soc_median_hourly_all), ONS_temp[rown, "median_all"],
                    soc_median_hourly_all) ]

        ##if even 2-digit info for soc_mean, soc_median missing, use soc_mean / median for all
        D[year==i,
            soc_mean_hourly:=ifelse(!is.na(soc) &
                is.na(soc_mean_hourly), soc_mean_hourly_all,
                    soc_mean_hourly)]
        D[year==i,
            soc_median_hourly:=ifelse(!is.na(soc) &
                is.na(soc_median_hourly), soc_median_hourly_all,
                    soc_median_hourly)]
        #convert to 2015 currency  (this syntax works for data.frame only, cpi)
        D[year==i, soc_mean_hourly:=soc_mean_hourly  / cpi[i-2001, 2] ]
        D[year==i, soc_median_hourly:=soc_median_hourly  / cpi[i-2001, 2] ]
}

D[, `:=`(
    soc_2d = factor(substr(SOC, 1, 2)), 
    year = factor(year, levels = c(2002:2016)), 
    age = factor(age, levels = c(35:64)) 
)]

X = model.matrix(~ ( age  + soc_2d + year ) *log(soc_mean_hourly)*log(soc_median_hourly)*male, data=D)

D[!is.na(SOC), log_y_hourly :=  X %*% coef] 
return(D)

}
