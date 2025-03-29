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

# Test if each are in the data, by testing the dictionary.
path = os.getcwd()
dictCsv = glob.glob(os.path.join(path, "*.data_dictionary.csv"))[0]
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
# Define aprent child by IBS_0 threshold, 0.0012 (above siblings, below parent).
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

# Get FID family identifier, by clustering on family networks.
import networkx as nx
famGraph = nx.Graph()
ID1List = relatedKinData["ID1"].tolist()
ID2List = relatedKinData["ID2"].tolist()
kinshipList = relatedKinData["Kinship"].tolist()
famEdges = [[ID1List[i], ID2List[i]] for i in range(0, len(relatedKinData))]
    #if kinshipList <= firstDegreeThreshold[1]]
# Add nodes to network (duplicates do not matter here), and their edge connections.
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
        # Define the family ID as the lowest ID of the family members.
        for famConnection in famConnections :
            relatedKinData.loc[((relatedKinData["ID1"] == famConnection) |
                (relatedKinData["ID2"] == famConnection)) & (
                    pd.isna(relatedKinData["FID"])), "FID"] = ID1

#!Test: are all FIDs filled in?
print(relatedKinData[pd.isna(relatedKinData["FID"])])
ID1 = np.random.choice(ID1List + ID2List)
print(ID1)
print(relatedKinData[(relatedKinData["ID1"] == ID1) | (relatedKinData["ID2"] == ID1)])

# Make family ID have no decimals
relatedKinData["FID"] = pd.to_numeric(
    relatedKinData["FID"], downcast="integer", errors="coerce")

# Get the parent/sibling identifier.
# InfType: Inferred relationship type, such as Dup/MZTwin, PO, FS, 2nd, 3rd, 4th, UN
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
    "kinship-adjusted.dat", sep=" ", index=False)


################################################################################
## Step 5. Create input files for Snipar imputation. 2. Age-sex file.
# https://snipar.readthedocs.io/en/latest/input%20files.html#agesex-file

# Get the base pheno data. file for Snipar imputation.
subprocess.check_call(["dx", "download", "data-clean/phenotype-extract.csv"])
phenoData = pd.read_csv("phenotype-extract.csv")
phenoKinData = phenoData.copy()

# A white-space delimited text file with header “FID”, “IID”, “sex”, “age”.
# Each row contains the family-ID, individual-ID, age, and sex of one individual.
# Male and Female sex should be represented with ‘M’ and ‘F’ respectively. 
# https://snipar.readthedocs.io/en/latest/input%20files.html#agesex-file
phenoKinData["IID"] = phenoKinData["participant.eid"]
phenoKinData["sex"] = ["M" if line == 1
    else "F"
    for line in phenoKinData["participant.p31"].tolist()]
phenoKinData["age"] = phenoKinData["participant.p21022"]

# Get the FAMID from the above clustering
phenoKinData["FID"] = phenoKinData["participant.eid"]


# Export for Snipar
phenoKinData[["FID", "IID", "sex", "age"]].to_csv(
    "phenotype-agesex.dat",
    index=False, sep=" ")
