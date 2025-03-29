#!/bin/bash
## Senan Hogan-Hennessy, 11 March 2025.
## Master bash script, for scripts showing IV without exclusion estimator.
# Note "R CMD BATCH" generates .Rout files showing consol output of each script.

################################################################################
## Scripts to be run on the UKB server, in "data/ukb-restricted/input/rep-code"
# Extract the relevant data files, from raw formats.
python3 data-extract.py
# Get the relatedness data
#Todo -> this is getting superceded by PRS_impute.
R CMD BATCH --no-save relatedness-interact.R
# Get the Ed PGI, using Okbay+ (2022) weights.
cd ../pgi-okbay-2022
R CMD BATCH --no-save tsv-convert.R
cd ..
bash prs-pipeline.sh
# Impute the parents Ed PGI, using sibling data (Snipar, Young+ 2022).
bash prs-impute.sh
# TODO


################################################################################
## Build data sets of relevance

# 1. HRS
cd hrs/data-build
# Combine panel of HRS data with genetic educ scores
R CMD BATCH --no-save hrs-build.R
# Collapse the panel of HRS data, keeping genetic educ scores
R CMD BATCH --no-save hrs-collapse.R
cd ../..

# 2. UKB
cd ukb/data-build
# Build the cross-section of phenotype data, with PGIs.
R CMD BATCH --no-save ukb-build.R
# Identify relatives in UKB data, based on genetic similarity.
#TODO -> this is becoming redundant with the Snipar imputing...
python3 -u ukb-related.py > ukb-related.log
# Go back to base folder.
cd ../..

## Statistical analysis.
cd data-analyse
# Summarise UKB data.
R CMD BATCH --no-save ukb-summarise.R
# Analsyse HRS data, with the invalid instrument educ score.
R CMD BATCH --no-save hrs-educ-iv.R





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
# cp paper.pdf ../direct-indirect-education-2025.pdf
