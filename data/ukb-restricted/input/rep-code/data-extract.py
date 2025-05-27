#!/usr/bin/python3
## Senan Hogan-Hennessy, 14 March 2025
## Extract columns from the phenotype data, given provided column names.
# https://dnanexus.gitbook.io/uk-biobank-rap/working-on-the-research-analysis-platform/accessing-data/accessing-phenotypic-data

# Packages
import dxpy
import subprocess
import os
import glob
# Install numpy, pandas, and networkx
subprocess.check_call(["pip", "install", "pandas", "networkx"])
import numpy as np
import pandas as pd
import networkx as nx


################################################################################
## Step 1. Get the file containing all available data-fields in phenotype data.

# Automatically discover dispensed dataset ID and load the dataset 
dispensed_dataset_id = dxpy.find_one_data_object(
    typename="Dataset", name="app*.dataset", folder="/", name_mode="glob")["id"]

# Get project ID
project_id = dxpy.find_one_project()["id"]

# Host raw dataset.
dataset = (":").join([project_id, dispensed_dataset_id])
print(dataset)
cmd = ["dx", "extract_dataset", dataset, "-ddd", "--delimiter", ","]
subprocess.check_call(cmd)


################################################################################
## Step 2. Collect all data-field names of interest, for pheno data.

# Input: list of variables I have parsed through, and want.
# See -> phenotype-names.txt
desiredFields = ["eid", # eid, person identifier (anonymised).
    "p31", # sex
    "p34", # birthyear
    "p52", # eid
    "p53_i0",   "p53_i1",     "p53_i2",    "p53_i3",   # Date visit (4 instances).
    "p738_i0",  "p738_i1",    "p738_i2",   "p738_i3",  # householdincome_cat (4 instances).
    "p21022",    # recruitedage
    "p20143",    # datelastcontact
    "p6138_i0",  "p6138_i1",  "p6138_i2",  "p6138_i3", # edquals (4 instances).
    "p845_i0",   "p845_i1",   "p845_i2",               # agefinishededuc (3 instances).
    "p21000_i0", # ethnicity
    "p22006",    # genetic_race
    "p20118_i0", # urbancat
    "p20277_i0", "p20277_i1", "p20277_i2", "p20277_i3", # jobcode (4 instances).
    "p6142_i0",  "p6142_i1",  "p6142_i2",  "p6142_i3", # employment (4 instances).
    "p767_i0",   "p767_i1",   "p767_i2",   "p767_i3",  # hours_workweek (4 instances). 
    "p816_i0",   "p816_i1",   "p816_i2",   "p816_i3",  # manualwork (4 instances). 
    "p129_i0",   "p129_i1",   "p129_i2",               # birth_n_coord (3 instances). 
    "p130_i0",   "p130_i1",   "p130_i2",               # birth_e_coord (3 instances). 
    "p1647_i0",  "p1647_i1",  "p1647_i2",              # birth_country (3 instances).
    "p40000_i0", # datedeath
    "p22020",    # inPCA
    "p26210",    # asthma_pgi
    "p26214",    # bipolar_pgi
    "p26216",    # bmi_pgi
    "p26240",    # height_pgi
    "p26275",    # schizophrenia_pgi
    "p26285"]    # t2diabetes_pgi

#TODO list, extracting new data:
# 1. Adjust the extract of employment category to be binary measures (not lists).
#    Same for urban rural designator.
# https://biobank.ctsu.ox.ac.uk/crystal/field.cgi?id=6142
# https://biobank.ctsu.ox.ac.uk/crystal/field.cgi?id=20118

# 3. Funky coding for birth coordinates -> Card uni proximity instrument.
# -> Get data to have yearly coordinates of UK higher ed (equivalent to IPEDS)
# -> Match these data to UKB, getting binary (in local area) and distance to closest
# -> gets columns for uni in local area at age 18, and distance to closest.

desiredFieldsStr = [f"participant.{f}" for f in desiredFields]
fieldNames = ",".join(desiredFieldsStr)
print(fieldNames)

# Test if each are in the data, by testing data dictionary.
dataPath = os.getcwd()
dictCsv = glob.glob(os.path.join(dataPath, "*.data_dictionary.csv"))[0]
dictData = pd.read_csv(dictCsv)
print([field + " not included." for field in desiredFields
    if field not in dictData["name"].tolist()])


################################################################################
## Step 3. Extract pheno data using dx extract_dataset command.

# Piece by piece data extraction (cannot get all in one call)
desiredLen = 9
desiredFieldsLists = [
    ["participant.eid"] + desiredFieldsStr[1:][i : i + desiredLen]
    for i in range(0, len(desiredFieldsStr[1:]), desiredLen)]
# Loop across each set of columns.
for i in range(0, len(desiredFieldsLists)) :
    print("Starting phenotype-extract-part" + str(i) + ".csv")
    desiredFieldNames = ",".join(desiredFieldsLists[i])
    subprocess.check_call(["dx", "extract_dataset",
        dataset, "--fields", desiredFieldNames, "--delimiter", ",", "--output",
        "phenotype-extract-part" + str(i) + ".csv"])
    print(i + 1, "out of", len(desiredFieldsLists),
        "datasets (column-wise) extracted.")

# Stitch back together the disparate datasets
phenoData = pd.read_csv("phenotype-extract-part0.csv")
for i in range(1, len(desiredFieldsLists)) :
    extraData = pd.read_csv("phenotype-extract-part" + str(i) + ".csv")
    phenoData = phenoData.merge(extraData, how="left", on="participant.eid")

## Clean the ed quals columns (needs python list comprehension, so do here).
# extract the column of lists, for each instance.
from ast import literal_eval
edQualsLists_i0 = [literal_eval(row)
    for row in phenoData["participant.p6138_i0"].fillna("[]").tolist()]
edQualsLists_i1 = [literal_eval(row)
    for row in phenoData["participant.p6138_i1"].fillna("[]").tolist()]
edQualsLists_i2 = [literal_eval(row)
    for row in phenoData["participant.p6138_i2"].fillna("[]").tolist()]
edQualsLists_i3 = [literal_eval(row)
    for row in phenoData["participant.p6138_i3"].fillna("[]").tolist()]
# Put these lists to one.
edQualsLists = [edQualsLists_i0[i] + edQualsLists_i1[i] +
    edQualsLists_i2[i] + edQualsLists_i3[i]
    for i in range(0, len(edQualsLists_i0))]
# Remove duplicates within each list 
edQualsLists = [list(dict.fromkeys(edQualsLists[i]))
    for i in range(0, len(edQualsLists_i0))]

# Get binary values for level of education.
# qual == -7    None of above, given mandatory min   -> edYears =  7    
# qual == -3  Prefer not to answer.                  -> edYears =  np.nan 
# qual ==  1  Uni degree.                            -> edYears =  20    
# qual ==  2  A Levels.                              -> edYears =  13    
# qual ==  3  CSEs (old version of exams age 16      -> edYears =  10    
# qual ==  4  GCSEs (new version of exams age 1      -> edYears =  10    
# qual ==  5  Vocational degree (one year higher ed) -> edYears =  19    
# qual ==  6  Professional qual (teach, nurse).      -> edYears =  15    
higheredList    = [int(1 in edQualsList) for edQualsList in edQualsLists]
alevelsList     = [int(2 in edQualsList) for edQualsList in edQualsLists]
gcsesList       = [int(3 in edQualsList or 4 in edQualsList)
    for edQualsList in edQualsLists]
vocationalList  = [int(5 in edQualsList) for edQualsList in edQualsLists]
professionaList = [int(6 in edQualsList) for edQualsList in edQualsLists]
edminimumList   = [int(-7 in edQualsList) for edQualsList in edQualsLists]
edmissingList   = [int(edQualsLists[i] == [] or edQualsLists[i] == [-3])
    for i in range(0, len(edQualsLists))]

# Put onto the phenotype dataframe.
phenoData["participant.edqual_highered"]     = higheredList
phenoData["participant.edqual_alevels"]      = alevelsList
phenoData["participant.edqual_gcses"]        = gcsesList
phenoData["participant.edqual_vocational"]   = vocationalList
phenoData["participant.edqual_professional"] = professionaList
phenoData["participant.edqual_minimum"]      = edminimumList
phenoData["participant.edqual_missing"]      = edmissingList

# Show that it worked.
print(phenoData[[
    "participant.eid",
    "participant.p6138_i0",
    "participant.edqual_highered",
    "participant.edqual_alevels",
    "participant.edqual_gcses",
    "participant.edqual_vocational",
    "participant.edqual_professional",
    "participant.edqual_minimum",
    "participant.edqual_missing"]])

# Save (with the adjusted columns).
phenoData.to_csv("phenotype-extract.csv")
# Upload to clean data output.
uploadCall = ["dx", "upload", "phenotype-extract.csv",
    "--path", "data-clean/phenotype-extract.csv"]
subprocess.check_call(uploadCall)

#! Test: this should have zero rows -> everyone has ed value, or marked as missing.
phenoData["participant.edqual_any"] = (
    phenoData["participant.edqual_highered"] +
    phenoData["participant.edqual_alevels"] +
    phenoData["participant.edqual_gcses"] +
    phenoData["participant.edqual_vocational"] +
    phenoData["participant.edqual_professional"] +
    phenoData["participant.edqual_minimum"])
phenoData[["participant.eid","participant.edqual_any"]][
    (phenoData["participant.edqual_any"] == 1) &
    (phenoData["participant.edqual_missing"] == 1)]


################################################################################
## Step 5. Extract relatedness file.

# Get the relatedness file.
# dx download "Bulk/Genotype Results/Genotype calls/ukb_rel.dat"
subprocess.check_call(["dx", "download",
    '"Bulk/Genotype Results/Genotype calls/ukb_rel.dat"'])

# Load the relatedness file.
relatedBaseData = pd.read_csv("ukb_rel.dat", sep = " ")

# Export TSV readable version.
relatedBaseData.to_csv("relatedness-extract.tsv",
    index=False, sep="\t")
subprocess.check_call(["dx", "upload", "relatedness-extract.tsv", "--path",
    '"data-clean/relatedness-extract.tsv"'])

# Show format of the file.
relatedData = relatedBaseData.copy()
print(relatedData)
## Compare the observation count to Muslimnova+ (2024) Table A1.
# Define 1st degree sibling/parents. (lower bound to take out identical twins).
firstDegreeThreshold = [0.1770, 0.3540]
relatedData["relative_1stdegree"] = (
    firstDegreeThreshold[0] < relatedData["Kinship"]) & (
        relatedData["Kinship"] <= firstDegreeThreshold[1])
# Define parent, child by IBS_0 threshold, 0.0012 (above siblings, below parent).
parentThreshold = 0.0012
relatedData["relative_parent"] = (
    relatedData["IBS0"] < parentThreshold)
# SHow the number of pairs in the sample.
print(relatedData[relatedData["relative_1stdegree"] == True].groupby(
    ["relative_1stdegree", "relative_parent"]).size())
# WANT, from Muslimnova+ (2024):
#                 |  1st degree /   |  1st degree
#                 |  Parent-child   |  siblings
# Kinship coef    |  0.1770- 0.3540 |  0.1770 -0.3540
# IBS 0           |  <0.0012        |  >0.0012
# N (pairs)       |  6,271          |  22,659


################################################################################
## Step 5. Create input files for Snipar imputation. 1. Kinship file
# https://snipar.readthedocs.io/en/latest/input%20files.html#kinship-file

# Get a copy of the base data.
relatedKinData = relatedData.copy()

# Get FID family identifier, by components on a family related network.
import networkx as nx
famGraph = nx.Graph()
ID1List = relatedKinData["ID1"].tolist()
ID2List = relatedKinData["ID2"].tolist()
kinshipList = relatedKinData["Kinship"].tolist()
famEdges = [[ID1List[i], ID2List[i]] for i in range(0, len(relatedKinData))]
    #if kinshipList <= firstDegreeThreshold[1]]
# Add nodes to network (duplicates do not matter here), + edge connections.
famGraph.add_nodes_from(ID1List)
famGraph.add_nodes_from(ID2List)
famGraph.add_edges_from(famEdges)
# The connected components are separate families.
relatedKinData["FID"] = np.nan
# Loop across the EIDs, filling in family EID, based on connected components.
i = 0
for ID1 in ID1List:
    print(str(i) + " out of " + str(len(ID1List)))
    i += 1
    famConnections = list(nx.node_connected_component(famGraph, ID1))
    if len(famConnections) == 1 :
        relatedKinData.loc[relatedKinData["ID1"] == ID1, "FID"] = ID1
    elif pd.notna(relatedKinData.loc[(relatedKinData["ID1"] == ID1), "FID"]).all() :
        continue
    elif len(famConnections) > 1 :
        # Define the family ID as the lowest ID among connected family.
        for famConnection in famConnections :
            relatedKinData.loc[((relatedKinData["ID1"] == famConnection) |
                (relatedKinData["ID2"] == famConnection)) & (
                    pd.isna(relatedKinData["FID"])), "FID"] = ID1

# Test: are all FIDs filled in?
print(relatedKinData[pd.isna(relatedKinData["FID"])])
ID1 = np.random.choice(ID1List + ID2List)
print(ID1)
print(relatedKinData[(relatedKinData["ID1"] == ID1) | (
    relatedKinData["ID2"] == ID1)])
# Make family ID have no decimals
relatedKinData["FID"] = pd.to_numeric(
    relatedKinData["FID"], downcast="integer", errors="coerce")

# Get the parent/sibling identifier.
# InfType: Inferred relationship type, such as Dup/MZTwin, PO, FS, 2nd, ...
relatedChoices = ["Dup/MZ", "PO", "FS"]
relatedConditions = [
    (firstDegreeThreshold[1] < relatedKinData["Kinship"]),
    (firstDegreeThreshold[0] < relatedKinData["Kinship"]) & (
        relatedKinData["Kinship"] <= firstDegreeThreshold[1]) & (
            relatedKinData["IBS0"] < parentThreshold),
    (firstDegreeThreshold[0] < relatedKinData["Kinship"]) & (
        relatedKinData["Kinship"] <= firstDegreeThreshold[1]) & (
            relatedKinData["IBS0"] >= parentThreshold)
]
relatedKinData["InfType"] = np.select(
    relatedConditions, relatedChoices, default="2nd")
# Show the number of different kind of relatives
from collections import Counter
print(Counter(relatedKinData["InfType"]))

# Export usable format, for PRS impute.
relatedKinData[["FID", "ID1", "ID2", "InfType"]].to_csv(
    "kinship-adjusted.dat", sep="\t", index=False)
# Upload to workspace
subprocess.check_call(["dx", "upload", "kinship-adjusted.dat", "--path",
    "data-input/kinship-adjusted.dat"])


################################################################################
## Step 5. Create input files for Snipar imputation. 2. Age-sex file.
# https://snipar.readthedocs.io/en/latest/input%20files.html#agesex-file

# Get the base pheno data. file for Snipar imputation.
subprocess.check_call(["dx", "download", "data-clean/phenotype-extract.csv"])
phenoData = pd.read_csv("phenotype-extract.csv")
phenoKinData = phenoData.copy()

# https://snipar.readthedocs.io/en/latest/input%20files.html#agesex-file
# A white-space delimited text file with header “FID”, “IID”, “sex”, “age”.
# Each row contains the family-ID, individual-ID, age, and sex of one individual.
# Male and Female sex should be represented with ‘M’ and ‘F’ respectively. 
phenoKinData["IID"] = phenoKinData["eid"]
phenoKinData["sex"] = ["M" if line == 1
    else "F"
    for line in phenoKinData["31-0.0"].tolist()]
phenoKinData["age"] = phenoKinData["21022-0.0"]

# Get the FAMID from the above clustering
famidData = pd.concat(
    [relatedKinData[["ID1", "FID"]].rename(columns={"ID1" : "IID"}),
        relatedKinData[["ID2", "FID"]].rename(columns={"ID2" : "IID"})],
            ignore_index=True).groupby(by=["IID", "FID"], as_index=False).size()
# Merge onto these data.
phenoKinData = phenoKinData.merge(famidData[["IID", "FID"]],
    how="left", on="IID")
# Missing values are people not even in the related data, because no relatives.
phenoKinData.loc[pd.isna(phenoKinData["FID"]), "FID"] = phenoKinData.loc[
    pd.isna(phenoKinData["FID"]), "IID"]
print(phenoKinData.loc[pd.isna(phenoKinData["FID"]), ])
# Make family ID have no decimals
phenoKinData["FID"] = pd.to_numeric(
    phenoKinData["FID"], downcast="integer", errors="coerce")

# Export for Snipar
phenoKinData[["FID", "IID", "sex", "age"]].to_csv(
    "phenotype-agesex.dat",
    index=False, sep="\t")
# Upload to workspace
subprocess.check_call(["dx", "upload", "phenotype-agesex.dat", "--path",
    "data-input/phenotype-agesex.dat"])
