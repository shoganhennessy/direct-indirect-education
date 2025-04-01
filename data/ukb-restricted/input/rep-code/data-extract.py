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
dataset = (":").join([project_id, dispensed_dataset_id])

cmd = ["dx", "extract_dataset", dataset, "-ddd", "--delimiter", ","]
subprocess.check_call(cmd)


################################################################################
## Step 2. Collect all data-field names of interest, for pheno data.

# Input: list of variables I have parsed through, and want.
# See -> phenotype-names.txt
desiredFields = ["eid",
    "p31",
    "p34",
    "p52",
    "p53_i0",
    "p738_i0",
    "p21022",
    "p20143",
    "p6138_i0",
    "p845_i0",
    "p21000_i0",
    "p22006",
    "p20277_i0",
    "p22599",
    "p22661",
    "p22601_a0",
    "p22617_a0",
    "p22602_a0",
    "p22604_a0",
    "p22605_a0",
    "p22501",
    "p40000_i0",
    "p22020"]
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

# Extract from the phenotype database.
cmd = ["dx", 
    "extract_dataset",
    dataset,
    "--fields",
    fieldNames,
    "--delimiter",
    ",",
    "--output",
    "phenotype-extract.csv"
]
subprocess.check_call(cmd)

os.chdir("/home/s/Dropbox/direct-indirect-education/data/ukb-restricted/input")

# Add on a column for Ed Years, to phenotype file.
phenoData = pd.read_csv("phenotype-extract.csv")

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
    for row in phenoData["participant.p6138_i0"].fillna("[]").tolist()]
# Get years of education for each qualification they have obtained.
edYearsLists = [
    [edQualsConvert(edQual) for edQual in edQualsList]
    for edQualsList in edQualsLists]
# Get the maximum edyears.
maxEdList = [
    max(edList)
    if len(edList) > 0
    else np.nan
    for edList in edYearsLists]
# Put onto the phenotype dataframe.
phenoData["participant.edyears"] = maxEdList
# Get the qualification corresponding to max Ed Years.
maxEdqualList = [
    edQualsLists[i][np.argmax(max(edYearsLists[i]))]
    if len(edQualsLists[i]) > 0
    else np.nan
    for i in range(0, len(edQualsLists))]
# Put onto the phenotype dataframe.
phenoData["participant.p6138_i0"] = maxEdqualList

# SHow that it worked.abs
print(phenoData[
    ["participant.eid", "participant.p6138_i0", "participant.edyears"]])

# Save (with the adjusted columns).
phenoData.to_csv("phenotype-extract.csv")
# Upload to clean data output.
uploadCall = ["dx", "upload", "phenotype-extract.csv",
    "--path", "data-clean/phenotype-extract.csv"]
subprocess.check_call(uploadCall)


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
for ID1 in ID1List :
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
    '"data-clean/kinship-adjusted.dat"'])


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
phenoKinData["IID"] = phenoKinData["participant.eid"]
phenoKinData["sex"] = ["M" if line == 1
    else "F"
    for line in phenoKinData["participant.p31"].tolist()]
phenoKinData["age"] = phenoKinData["participant.p21022"]

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
    '"data-clean/phenotype-agesex.dat"'])
