#!/usr/bin/R
## Senan Hogan-Hennessy, 28 Feb 2025
## Extract the UKB relatedness file, on UKB platform.

# Get the relatedness file.
source('dx download "Bulk/Genotype Results/Genotype calls/ukb_rel.dat"')
# Show what we are dealing with.
print(readLines("ukb_rel.dat", n = 10))

# Load the relatedness file.
relatedness.data <- read.table("ukb_rel.dat", header = TRUE)

# Show format of the file.
print(relatedness.data[1:10, ])

# Export usable data file.
write.table(relatedness.data, "relatedness-extract.tsv",
    row.names = FALSE, sep = "\t")
source(paste('dx upload relatedness-extract.tsv --path',
    '"data-clean/relatedness-extract.tsv"'))
