#!~/venv/bin/python3
## Senan Hogan-Hennessy, 18 March 2025
## Define UKB parents + siblings, from provided relatedness file.

## Packages:
# Package for data files.
import pandas as pd
# Package for numeric manipulation
import numpy as np
# Package for file system navigation.
import os

## Set up environment.
#! Temporary ->
#! Change to file location (VSCode makes this impossible to automate locally...)
os.chdir(os.path.join(
    "/home", "s", "Dropbox", "direct-indirect-education", "programs", "ukb", "data-build"))
# Local NP seed set
rng = np.random.default_rng(47)
# Show operational time, date.
import datetime
time = datetime.datetime.now()
print(time.strftime("%H:%M:%S, %d %B %Y"))
# Define local of input files.
dataFolder = os.path.join("..", "..", "..", "data", "ukb-restricted")
inputFolder = os.path.join(dataFolder, "input")
intermedFolder = os.path.join(dataFolder, "intermediate-files")


################################################################################
## Load UKB data.

# Load the pre-cleaned phenotype date, for DOB and gender.
phenoData = pd.read_csv(
    os.path.join(intermedFolder, "ukb-intermed-pheno.csv"),
    sep=",", header=0)
# Set the index as EID
phenoData.set_index(phenoData["eid"], drop=False, inplace=True)

# Load the raw relatedness file.
relatedBaseData = pd.read_csv(
    os.path.join(inputFolder, "relatedness-extract.tsv"),
    sep="\t", header=0)


################################################################################
## Clean and extract the phenotype  education data (Python lists comprehension).

# extract the column of lists.
from ast import literal_eval
edqualsLists = [literal_eval(row)
    for row in phenoData["edquals"].fillna("[]").tolist()]
# Get the maximum qualification.
maxEdList = [max(edList)
    if len(edList) > 0
    else np.nan
    for edList in edqualsLists]
# Convert ed quals -> edyears, Following ISCED (see Mulsimnova+ 2024, p14)
# https://github.com/DilnozaM/Rank-Concordance-of-PGI/blob/main/GWAS%26PGS/EA_GWAS/Preparing_residual_EA_new_nosibrel_fastgwas.do
edyearsList = [
    np.nan      if maxEdQual == np.nan
    else 7      if maxEdQual == -7 # None of above, so code mandatory minimum.
    else np.nan if maxEdQual == -3 # Prefer not to answer.
    else 20     if maxEdQual == 1  # Uni degree.
    else 13     if maxEdQual == 2  # A Levels.
    else 10     if maxEdQual == 3  # CSEs (old version of exams at age 16)
    else 10     if maxEdQual == 4  # GCSEs (new version of exams at age 16)
    else 19     if maxEdQual == 5  # Vocational degree (one year of tertiary).
    else 15     if maxEdQual == 6  # Professional qualification (teach, nurse).
    else np.nan
    for maxEdQual in maxEdList]
# Put onto the phenotype dataframe.
phenoData["edyears"] = edyearsList


################################################################################
## Identify relatives from UKB relatedness file.

# King manual threshold for 1st degree related (sibling or parent), >= 0.1770.
relatedData = relatedBaseData[relatedBaseData["Kinship"] >= 0.1770]

# Sort within each person (ID1) in descending relatedness, and set index
relatedData = relatedData.sort_values(
    by=["ID1", "Kinship"], ascending=[True, False])
relatedData["eid"] = relatedData["ID1"]
relatedData.set_index(relatedData["eid"], drop=True, inplace=True)

# For each individual (ID1), get the list of 1st degree relatives (ID2).
# Note: list of relatives is in descending relatedness order, i.e. [1st, 2nd, .]
relatedDict = {}
for id1 in relatedData["ID1"].tolist() :
    if id1 in relatedDict :
        continue
    else :
        relatedDict[id1] = relatedData[
            relatedData["ID1"] == id1]["ID2"].tolist()

# Get an empty data frame, to fill, and return.
parentData = pd.DataFrame({}, columns=["eid", "eid_mother", "eid_father",
    "eid_sibling1", "eid_sibling2", "eid_sibling3"])
parentData["eid"] = list(relatedDict.keys())
parentData.set_index(parentData["eid"], drop=False, inplace=True)

# Identify each relative, by gender + age difference, and add to parentData.
for id1 in list(relatedDict.keys()) :
    ownData = phenoData.loc[id1]
    for relativeEid in relatedDict[id1] :
        relativeData = phenoData.loc[relativeEid]
        # Test if they could be child-parent, 15 <= age gap <= 60
        ageGap = ownData["birthyear"] - relativeData["birthyear"]
        if 15 <= ageGap & ageGap <= 60 :
            # Parent, who is male, is father.
            if relativeData["sex_male"] == 1 :
                parentData.loc[id1, "eid_father"] = relativeEid
            # Parent, who is female, is mother.
            elif relativeData["sex_male"] == 0 :
                parentData.loc[id1, "eid_mother"] = relativeEid
            # Then exit this loop
            continue
        # If not child-parent, then siblings
        elif ageGap <= 15 :
            if pd.isna(parentData.loc[id1, "eid_sibling1"]) :
                parentData.loc[id1, "eid_sibling1"] = relativeEid
            elif pd.isna(parentData.loc[id1, "eid_sibling2"]) :
                parentData.loc[id1, "eid_sibling2"] = relativeEid
            elif pd.isna(parentData.loc[id1, "eid_sibling3"]) :
                parentData.loc[id1, "eid_sibling3"] = relativeEid

# Show how many people have any 1st degree relatives -> all in the Kinship file.
print(parentData[pd.notna(parentData["eid_father"])
    | pd.notna(parentData["eid_mother"])
    | pd.notna(parentData["eid_sibling1"])])


################################################################################
## Merge relative EIDs onto the phenotype file.

#! Todo:
#! Also figure out the imputation steps.
