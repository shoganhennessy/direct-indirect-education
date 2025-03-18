#!/bin/bash
## Senan Hogan-Hennessy, 11 March 2025.
## Master bash script, for scripts showing IV without exclusion estimator.
# Note "R CMD BATCH" generates .Rout files showing consol output of each script.

## Build data sets of relevance
cd data-build

# 1. HRS
# Combine panel of HRS data with genetic educ scores
R CMD BATCH --no-save hrs-build.R
# Collapse the panel of HRS data, keeping genetic educ scores
R CMD BATCH --no-save hrs-collapse.R

# 2. UKB
# Build the cross-section of phenotype data, with PGIs.
R CMD BATCH --no-save ukb-build.R
# Identify relatives in UKB data, based on genetic similarity.
python3 ukb-related.py
# Go back to base folder.
cd ..

## Statistical analysis.
cd data-analyse
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
