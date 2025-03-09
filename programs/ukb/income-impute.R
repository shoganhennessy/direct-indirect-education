#!/usr/bin/R
## Senan Hogan-Hennessy, 18 Feb 2025
## Convert the Kweon+ (2024) UKB income imputation file to a usable format.
# Source: 1. https://doi.org/10.62891/aac85602
#         2. https://doi.org/10.1038/s41562-024-02080-7

## Set up the R environment
set.seed(47)
# Show the date:
print(format(Sys.time(), "%H:%M %Z %A, %d %B %Y"))

## Load libraries
# Functions for data manipulation and visualisation
library(tidyverse)
# Define where data files are.
data.folder <- file.path("..", "..", "data", "earnings-imputed")


################################################################################
## Load and clean the Kweon+ (2015) provided data file.
load(file.path(data.folder, "data_input.Rdata"))

# Gives a list of data files ``ONSlist'' -> harmonise into one data file.
#! Note: these wages are nominal hourly figures, for years 2002--2016.
#! Note: and they need updating to 2025 pounds (British CPI).
