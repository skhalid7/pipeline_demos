#!/usr/bin/env python3

import pandas as pd
import numpy as np
import sys, re
from collections import defaultdict

def makeInfoCol(input_cols, col_names):
    return_strings = []
    for col in col_names:
        return_strings.append("=".join([col, str(input_cols[col])]))
    return ";".join(return_strings)

def ACGTonly(allele):
    return len([ x for x in allele if x not in ["A", "C", "G", "T"] ]) == 0

#GLOBAL VARIABLES
COLS_2_KEEP_FROM_INPUT = ["CHR", "POS", "ID", "otherAllele", "effectAllele", "SE", "beta", "pvalue", "infoScore", "effect_AF", "cases_AF", "controls_AF", "snp_N"]
STATIC_FIELDS = ["total_N", "cases_N", "controls_N"]
COLS_2_KEEP_IN_INFO_COL = COLS_2_KEEP_FROM_INPUT[4:len(COLS_2_KEEP_FROM_INPUT)] + STATIC_FIELDS + ["phenotype", "study"]
DELIM = "\t"

#INPUT PARAMS
file_name = sys.argv[1]
output_base_path = sys.argv[2]

#Read in DD from Google Sheets
google_sheet_id = '1H2ePJPADlkgw1fsmHZaL_9BqPtl1U_LT1OLi8TSHf7Q'
sheet_name = "External_Summ_Stats"
url_sheet = r'https://docs.google.com/spreadsheets/d/{0}/gviz/tq?tqx=out:csv&sheet={1}'.format(google_sheet_id, sheet_name)

#Format DD, create array and dict to pick and rename relevant columns
df_data_dictionary  = pd.read_csv(url_sheet, sep = ",").query("phenotype_file == '{0}'".format(file_name))
if df_data_dictionary.shape[0] == 0:
    raise ValueError('Unable to find specified file in path, check file name and location')
if df_data_dictionary.shape[0] > 1:
    raise ValueError('Multiple files with the same name exist, check meta file')
df_data_dictionary.to_csv(output_base_path + "_params_file.tsv", sep = "\t", index = False)
col_numeric_indices = [i for i in range(0,len(df_data_dictionary.columns)) if df_data_dictionary.columns[i] in COLS_2_KEEP_FROM_INPUT]
rename_dict = df_data_dictionary.fillna("NA").iloc[0,col_numeric_indices].to_dict()
rename_dict = {v: k for k, v in rename_dict.items() if v != "NA"}
cols_to_pick = [key for( key, values) in rename_dict.items() if key != "NA"]
static_field_values = df_data_dictionary[STATIC_FIELDS].values[0]

#Read in input file and format to harmonized format
input_df = pd.read_csv(file_name, sep = DELIM, usecols = cols_to_pick,  na_values= [ -99, -999 ])
input_df.rename(columns = rename_dict, inplace=True)
for i in range(0,len(STATIC_FIELDS)):
    input_df[STATIC_FIELDS[i]] = static_field_values[i]
input_df["study"] = df_data_dictionary.study.values[0]
input_df["phenotype"] = df_data_dictionary.phenotype.values[0]

'''
QC Summary statistics 
'''

print(str(input_df.shape[0]) + " rows read in")
#effect and other allele columns to upper
input_df["effectAllele"] = input_df.effectAllele.apply(lambda x : x.upper())
input_df["otherAllele"] = input_df.otherAllele.apply(lambda x : x.upper())

input_df["ea_is_acgt"] = input_df.effectAllele.apply(lambda x : ACGTonly(x))
input_df["oa_is_acgt"] = input_df.otherAllele.apply(lambda x : ACGTonly(x))

#effect and other allele columns should be proper variants:
input_df = input_df.query("oa_is_acgt & ea_is_acgt")

#remove variants with null betas and SEs
input_df = input_df[input_df.beta.notnull() & input_df.SE.notnull()]

if "effect_AF" in input_df:
    #filter out variants with maf of 0.01%
    input_df = input_df.query("effect_AF >= 0.001 & effect_AF <= 0.999")
if "cases_AF" in input_df:
    #filter out variants where the cases maf is less than 0.01%
    input_df.query("cases_AF >= 0.001 & cases_AF <= 0.999")
if "cases_AF" in input_df and "cases_N" in input_df:
    #filter out variants where the expectes Allele Count < 25
    input_df["AC"] = input_df["cases_AF"] * input_df["cases_N"] * 2
    input_df = input_df.query("AC > 25")
if "infoScore" in input_df.columns:
    #filter out variants with info score <= 0.3
    input_df = input_df.query("infoScore > 0.3")

print(str(input_df.shape[0]) + " rows being formatted and written out")

#Create Info Col and other columns for VCF
info_cols = [col for col in COLS_2_KEEP_IN_INFO_COL if col in input_df.columns]
input_df["INFO"] = input_df[info_cols].apply(lambda x : makeInfoCol(x, info_cols), axis = 1)
input_df["QUAL"] = "."
input_df["FILTER"] = "."
input_df["POS"] = input_df.POS.astype(int)
input_df["ID"] = input_df["ID"].fillna(".") #fill missing id columns with '.'

#write out body of VCF file
input_df[["CHR", "POS", "ID", "otherAllele", "effectAllele", "FILTER", "QUAL", "INFO"]].to_csv("vcf_formatted_stats", sep = "\t", index = False, header = False)
