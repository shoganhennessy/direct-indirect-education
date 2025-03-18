#!/bin/python3
## Senan Hogan-Hennessy, 14 March 2025
## Extract columns from the phenotype data, given provided column names.
# https://dnanexus.gitbook.io/uk-biobank-rap/working-on-the-research-analysis-platform/accessing-data/accessing-phenotypic-data

# Packages
import dxpy
import subprocess
import os
import glob
# pip install -I pandas==1.3.5
import pandas as pd


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
## Step 2. Collect all data-field names of interest. 

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
## Step 3. Extract the data using dx extract_dataset command.

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
