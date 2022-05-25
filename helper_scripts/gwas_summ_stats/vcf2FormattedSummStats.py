#!/usr/bin/env python3

import pandas as pd
import numpy as np
import sys, re
from decimal import Decimal

def infoCol2Dict(info_col, cols_to_break):
    info_col = info_col.split(";")
    tuples = [(col.split("=")[0],col.split("=")[1]) for col in info_col if len(col.split("=")) == 2 and col.split("=")[0] in cols_to_break]
    return dict(tuples)

input_file = sys.argv[1]
lines_to_skip = int(sys.argv[2])
output_file = input_file.split(".vcf")[0] + ".summary_stats.formatted.tsv"
COLS_IN_FORMATTED_FILE = ["CHR", "POS", "ID", "REF", "ALT", "effectAllele", "OriginalContig", "OriginalStart", "SwappedAlleles", "effectAllele_pre_swap", "beta", "SE", "pvalue", "infoScore", "effect_AF", "cases_AF", "controls_AF", "total_N", "cases_N", "controls_N", "snp_N", "phenotype", "study"]

input_vcf = pd.read_csv(input_file, sep = "\t", skiprows = lines_to_skip)
input_vcf.rename(columns = {"#CHROM" : "CHR"}, inplace=True)
#filter out alt_contigs
input_vcf["CHR_check"] = input_vcf["CHR"].apply(lambda x : len(str(x).split("_")) == 1)
input_vcf = input_vcf.query("CHR_check")

#convert INFO field to dict
input_vcf["INFO_DICT"] = input_vcf["INFO"].apply(lambda x : infoCol2Dict(x, COLS_IN_FORMATTED_FILE))
#create each field in the dict into a new column
cols_to_make_from_df = input_vcf["INFO_DICT"].apply(lambda x : list(x.keys())).iloc[0:1].values[0]
for col in cols_to_make_from_df:
    input_vcf[col] = input_vcf["INFO_DICT"].apply(lambda x : x[col])
#for fields not in INFO add them in as missing
missing_cols_to_add = [col for col in COLS_IN_FORMATTED_FILE if col not in cols_to_make_from_df and col not in input_vcf.columns]
for col in missing_cols_to_add:
    input_vcf[col] = None

#add column for swapped alleles
input_vcf["SwappedAlleles"] = input_vcf["INFO"].apply(lambda x : "SwappedAlleles" in x)
input_vcf = input_vcf[COLS_IN_FORMATTED_FILE]
#Change beta, effect_AF, cases_AF, controls_AF and effectAllele column for swapped Alleles
subset_df = input_vcf.query("SwappedAlleles")
subset_df["beta"] = subset_df["beta"].apply(lambda x : Decimal(x) * -1)
for af_col in [col for col in ["effect_AF", "cases_AF", "controls_AF"] if col in cols_to_make_from_df]:
    subset_df[af_col] = subset_df[af_col].apply(lambda x : 1 - Decimal(x))
subset_df["effectAllele_pre_swap"] = subset_df["effectAllele"]
subset_df["effectAllele"] = subset_df["ALT"]
input_vcf.update(subset_df)

#type case POS to int (gets read in as double)
input_vcf["CHR"] = input_vcf["CHR"].apply(lambda x : str(x).split(r'.')[0])
input_vcf["POS"] = input_vcf.POS.astype(int)
#Create uid column
input_vcf["uid"] = input_vcf[["CHR", "POS", "REF", "ALT"]].apply(lambda cols : ":".join([str(x) for x in cols]), axis = 1)
input_vcf = input_vcf[["uid"] + COLS_IN_FORMATTED_FILE]
#filter out cases where contig changed between builds and new contig is an X chromosome contig
input_vcf = input_vcf.query("OriginalContig == CHR or CHR != 'X'")
#filter out multi-allelic / lift over sending two variants to the same position error
input_vcf = input_vcf.drop_duplicates(subset = 'uid', keep = False) #change this to locus, not to uid

#write out file
input_vcf.to_csv(output_file, sep = "\t", index = False)
