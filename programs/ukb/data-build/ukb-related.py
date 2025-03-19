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

# Define function edQual -> Edyears Following ISCED (see Mulsimnova+ 2024, p14)
# https://github.com/DilnozaM/Rank-Concordance-of-PGI/blob/main/GWAS%26PGS/EA_GWAS/Preparing_residual_EA_new_nosibrel_fastgwas.do
def edQualsConvert(edQual) :
    """ Input: edQual, a number representing ed qual from UKB
    Output: edYears, a number representing num years of ed from that qual.
    """
    edYears = np.nan
    if edQual == -7   : edYears =  7      # None of above, given mandatory min.
    elif edQual == -3 : edYears =  np.nan # Prefer not to answer.
    elif edQual == 1  : edYears =  20     # Uni degree.
    elif edQual == 2  : edYears =  13     # A Levels.
    elif edQual == 3  : edYears =  10     # CSEs (old version of exams age 16)
    elif edQual == 4  : edYears =  10     # GCSEs (new version of exams age 16)
    elif edQual == 5  : edYears =  19     # Vocational degree (one year higher).
    elif edQual == 6  : edYears =  15     # Professional qual (teach, nurse).
    return(edYears)
# extract the column of lists.
from ast import literal_eval
edQualsLists = [literal_eval(row)
    for row in phenoData["edquals"].fillna("[]").tolist()]
# Get years of education for each qualification they have obtained.
edYearsLists = [
    [edQualsConvert(edQual) for edQual in edQualsList]
    for edQualsList in edQualsLists]
# Get the maximum edyears.
maxEdList = [max(edList)
    if len(edList) > 0
    else np.nan
    for edList in edqualsLists]
# Put onto the phenotype dataframe.
phenoData["edyears"] = maxEdList


################################################################################
## Identify relatives from UKB relatedness file.

# King manual threshold for 1st degree related (sibling or parent), >= 0.1770.
print(relatedBaseData)
firstDegreeThreshold = [0.1770, 0.3540]
relatedData = relatedBaseData[
    (firstDegreeThreshold[0] < relatedBaseData["Kinship"]) & 
    (relatedBaseData["Kinship"] <= firstDegreeThreshold[1])]

# Sort within each person (ID1) in descending relatedness, and set index
relatedData = relatedData.sort_values(
    by=["ID1", "Kinship"], ascending=[True, False])
relatedData["eid"] = relatedData["ID1"]
relatedData.set_index(relatedData["eid"], drop=True, inplace=True)

# For each individual (ID1), get the list of 1st degree relatives (ID2).
# Note: list of relatives is in descending relatedness order, i.e. [1st, 2nd, .]
relatedDict = {}
for id1 in relatedData["ID1"].tolist() :
    if id1 not in relatedDict :
        relatedDict[id1] = relatedData[
            relatedData["ID1"] == id1]["ID2"].tolist()

# Get an empty data frame, to fill, and return.
parentData = pd.DataFrame({}, columns=["eid", "eid_mother", "eid_father",
    "eid_sibling1", "eid_sibling2", "eid_sibling3"])
parentData["eid"] = list(relatedDict.keys())
parentData.set_index(parentData["eid"], drop=False, inplace=True)

# Identify each relative, by gender + age difference, and add to parentData.
for id1 in parentData["eid"].tolist() :
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
        # If not child-parent, then siblings
        else :
            if pd.isna(parentData.loc[id1, "eid_sibling1"]) :
                parentData.loc[id1, "eid_sibling1"] = relativeEid
            elif pd.isna(parentData.loc[id1, "eid_sibling2"]) :
                parentData.loc[id1, "eid_sibling2"] = relativeEid
            elif pd.isna(parentData.loc[id1, "eid_sibling3"]) :
                parentData.loc[id1, "eid_sibling3"] = relativeEid

# Show how many people have any 1st degree relatives -> all in the Kinship file.
print(parentData)
print(parentData[pd.notna(parentData["eid_father"])
    | pd.notna(parentData["eid_mother"])
    | pd.notna(parentData["eid_sibling1"])])
# Show how many people have any parent + sibling degree relatives.
print(parentData[pd.notna(parentData["eid_father"])
    | pd.notna(parentData["eid_mother"])])
print(parentData[(pd.notna(parentData["eid_father"])
    | pd.notna(parentData["eid_mother"]))
    & pd.notna(parentData["eid_sibling1"])])

#TODO: Muslimnova (2024) uses IBS_0 < 0.0012 to define first-degree parents.
#TODO  and in so doing get a higher obs count. 
#TODO  Carvalho (2024) get higher match by another approach too.
#TODO  -> Needs refining, to get the most amount of parents.


################################################################################
## Merge relative EIDs onto the phenotype file.

# Define variables to take over
relativeVars = ["eid", "edpgi_norm"]
fatherVars   = {var : var + "_father"   for var in relativeVars}
motherVars   = {var : var + "_mother"   for var in relativeVars}
sibling1Vars = {var : var + "_sibling1" for var in relativeVars}
sibling2Vars = {var : var + "_sibling2" for var in relativeVars}
sibling3Vars = {var : var + "_sibling3" for var in relativeVars}

# Merge father Ed PGI onto the parentData.
edpgiData = phenoData[relativeVars].reset_index(drop=True).rename(
    columns=fatherVars)
parentData = parentData.merge(edpgiData, how="left", on="eid_father")
# Merge mother Ed PGI onto the parentData.
edpgiData = phenoData[relativeVars].reset_index(drop=True).rename(
    columns=motherVars)
parentData = parentData.merge(edpgiData, how="left", on="eid_mother")
# Merge sibling1 Ed PGI onto the parentData.
edpgiData = phenoData[relativeVars].reset_index(drop=True).rename(
    columns=sibling1Vars)
parentData = parentData.merge(edpgiData, how="left", on="eid_sibling1")
# Merge sibling2 Ed PGI onto the parentData.
edpgiData = phenoData[relativeVars].reset_index(drop=True).rename(
    columns=sibling2Vars)
parentData = parentData.merge(edpgiData, how="left", on="eid_sibling2")
# Merge sibling3 Ed PGI onto the parentData.
edpgiData = phenoData[relativeVars].reset_index(drop=True).rename(
    columns=sibling3Vars)
parentData = parentData.merge(edpgiData, how="left", on="eid_sibling3")

# Add on the relative variables to the phenotype data.
phenoData = phenoData.reset_index(drop=True).merge(
    parentData, how="left", on="eid")


################################################################################
## Save resulting data file.


df.to_csv(index=False)