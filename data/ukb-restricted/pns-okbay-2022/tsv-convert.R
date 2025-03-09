#!/usr/bin/R
## Senan Hogan-Hennessy, 28 Feb 2025
## Put the TSV GWAS summary stats file into CSV, with same columns as needed.
# Set up the R environment
set.seed(47)
# Show the date:
print(format(Sys.time(), "%H:%M %Z %A, %d %B %Y"))
# Load functions for data manipulation.
library(tidyverse)

# Read the provided TSV file, from SSGAC.
tsv.data <- read.table("EA4_additive_p1e-5_clumped.txt", header = TRUE)

# Get a version that mimics the format needed for the PRS pipeline.
csv.data <- tsv.data %>%
    transmute(
        rsid = rsID,
        chr_name = Chr,
        chr_position = BP,
        effect_allele = Effect_allele,
        noneffect_allele = Other_allele,
        Beta = Beta,
        eaf = EAF_HRC) %>%
    mutate(chr_pos = paste(chr_name, chr_position, sep = ":"))

# Write to a CSV file.
csv.data %>% write_csv("EA4_additive_p1e-5_clumped.csv")
