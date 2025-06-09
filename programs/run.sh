#!/bin/bash
## Senan Hogan-Hennessy, 11 March 2025.
## Master bash script, for scripts showing IV without exclusion estimator.
# Note "R CMD BATCH" generates .Rout files showing consol output of each script.

################################################################################
## Scripts to be run on the UKB server, in "data/ukb-restricted/input/rep-code"

# Extract the relevant data files, from raw formats.
python3 -u data-extract.py > data-extract.log

# Calculate the Ed PGI, using Okbay+ (2022) weights and prs-pipeline.
mkdir pgi-okbay-2022
cd pgi-okbay-2022
R CMD BATCH --no-save tsv-convert.R
bash prs-pipeline.sh "EA4_additive_p1e-5_clumped.csv" "raw-ed-pgi-okbay-2022"
bash prs-impute.sh "EA4_additive_p1e-5_clumped.csv" "imputed-ed-pgi-okbay-2022"
cd ..

# Calculate the Okbay+ (2022) PGI, weights from the excluded subsample.
cd pgi-okbay-excluded-2022
R CMD BATCH --no-save tsv-convert.R
bash prs-pipeline.sh "EA4_additive_excl_23andMe.csv" "raw-ed-pgi-okbay-exclude-2022"
bash prs-impute.sh "EA4_additive_excl_23andMe.csv" "imputed-raw-ed-pgi-okbay-exclude-2022"
cd ..

# Calculate the separate Ed PGI, using Tan+ (2024) weights and prs-pipeline.
mkdir pgi-tan-2024
cd pgi-tan-2024
dx download data-input/tan_2024_educational_attainment.tsv
# R CMD BATCH --no-save tsv-convert.R
dx download data-input/tan_2024_educational_attainment.csv
bash prs-pipeline.sh "tan_2024_educational_attainment.csv" "raw-ed-pgi-tan-2024"
bash prs-impute.sh "tan_2024_educational_attainment.csv" "imputed-ed-pgi-tan-2024"
cd ..


################################################################################
## Build data sets of relevance

cd ukb/data-build
# Build the cross-section of phenotype data, with PGIs (raw + parent imputed).
R CMD BATCH --no-save ukb-build.R
# Go back to base folder.
cd ../..

## Statistical analysis, UKB
cd ukb/data-analyse
# Summarise UKB data.
R CMD BATCH --no-save ukb-summarise.R
# First-stage effects, justifying the imputation research design.
R CMD BATCH --no-save ukb-first-stage.R
# Total genetic effect estimates.
R CMD BATCH --no-save ukb-genetic-effects.R
# Mediation analysis, direct + indirect effects of Ed PGI.
R CMD BATCH --no-save ukb-direct-indirect.R
# Go back to base folder.
cd ../..

## Statistical analysis, HRS
cd hrs/data-analyse

# Summarise HRS data.
R CMD BATCH --no-save hrs-summarise.R
# Analsyse HRS data, with the invalid instrument educ score.
R CMD BATCH --no-save hrs-educ-iv.R
# Go back to base folder.
cd ..

# Compile paper, with outputs from R scripts above as inputs for TeX files
cd ../text
latexmk -pdf paper.tex
latexmk -c
cp paper.pdf ../shoganhennessy-jmp-2025.pdf
