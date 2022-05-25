#!/usr/bin/env python

import pyspark.sql.functions as f
import pandas as pd
from pyspark.sql.types import StringType, BooleanType

'''
Set of functions to wrangle pvcfs using Glow
'''

def uid2GlowFmt(df):
    '''
    Converts a dataframe with C:P:R:A to a Glow format VCF for writing out or running through the pipe transformer
    '''
    return \
    df.withColumn("uid_col", f.split(f.col("uid"), ":"))\
    .withColumn("contigName", f.col("uid_col")[0])\
    .withColumn("start", (f.col("uid_col")[1] - 1).cast("Long"))\
    .withColumn("end", (f.col("start") + f.length(f.col("uid_col")[2])).cast("Long"))\
    .withColumn("referenceAllele", f.col("uid_col")[2])\
    .withColumn("alternateAlleles", f.array(f.col("uid_col")[3]))\
    .drop("uid_col")

def glowFmt2uid(df):
    '''
    Converts Glow VCF format to C:P:R:A format
    IMPORTANT ASSUMES SITE IS BIALLELIC, IF NOT run split_multi first
    '''
    return \
    df.withColumn("CHR", f.col("contigName"))\
    .withColumn("POS", f.col("start") + 1)\
    .withColumn("REF", f.col("referenceAllele"))\
    .withColumn("ALT", f.col("alternateAlleles")[0])\
    .withColumn("uid", f.concat_ws(":", f.array(f.col("CHR"), f.col("POS"), f.col("REF"), f.col("ALT"))))


def processGenotypes(df, extra_fields = []):
    '''
    1) Given a list of fields from the genotypes column to keep (default: calls, phased, depth, conditionalQuality, alleleDepths) removes extra fields from the genotype column (extra fields with array value can be problematic downstream)
    2) Recodes contigNames form chr12/CHR12/Chr12 -> 12
    3) creates an array "called" to keep track of missing calls, this array is used to fix multiallelics after being split.
    4) drops attributes column (unflattened info column)
    5) convert [1,-1] genotypes to [-1,-1] otherwise they get counted as homAA variants
    '''
    fields_to_keep = ["calls" ,"phased", "depth", "conditionalQuality", "alleleDepths"]
    fields_to_keep.extend(extra_fields)
    
    return df.withColumn("contigName", f.when(f.col("contigName").contains("chr"), f.split(f.col("contigName"), "r")[1])\
    .when(f.col("contigName").contains("CHR"), f.split(f.col("contigName"), "R")[1])\
    .when(f.col("contigName").contains("Chr"), f.split(f.col("contigName"), "r")[1])\
    .otherwise(f.col("contigName")))\
    .withColumn("genotypes", f.expr("transform(genotypes, g -> struct(if(array_contains(g.calls, -1), array(-1,-1), g.calls) as calls, g.sampleId as sampleId, g.phased as phased, g.depth as depth, g.conditionalQuality as conditionalQuality, g.alleleDepths as alleleDepths))"))\
    .withColumn("called", f.expr("transform(genotypes, g -> ! array_contains(g.calls, -1))"))\
    .drop("attributes")

def correctMultiAllelics(df):
    '''
    This function makes the split multiallelic command mimic hails split_multi function
    Using the "called" column created earlier will set missing calls in multi allelic sites to homRR if the other allele is called.
    HetAA (1/2) sites are converted to hetRA (0/1) for each allele
    '''

    return df.withColumn("genotypes_called", f.when(f.col("splitFromMultiAllelic"), f.arrays_zip(f.col("genotypes"), f.col("called"))))\
    .withColumn("genotypes_called", f.when(f.col("splitFromMultiAllelic"), f.expr("transform(genotypes_called, g -> if(g.called AND array_contains(g.genotypes.calls, 0) AND array_contains(g.genotypes.calls, -1), array(0,0),\
if(array_contains(g.genotypes.calls, 1) AND array_contains(g.genotypes.calls, -1) AND g.called, array(0,1), g.genotypes.calls)))")))\
    .withColumn("genotypes_called", f.when(f.col("splitFromMultiAllelic"), f.arrays_zip(f.col("genotypes"), f.col("genotypes_called"))))\
    .withColumn("genotypes", f.when(f.col("splitFromMultiAllelic"), \
    f.expr("transform(genotypes_called, g -> struct(g.genotypes_called as calls, g.genotypes.sampleId as sampleId, g.genotypes.phased as phased,g.genotypes.depth as depth, g.genotypes.conditionalQuality as conditionalQuality, g.genotypes.alleleDepths as alleleDepths))"))\
    .otherwise(f.expr("transform(genotypes, g -> struct(g.calls as calls, g.sampleId as sampleId, g.phased as phased,g.depth as depth, g.conditionalQuality as conditionalQuality, g.alleleDepths as alleleDepths))")))\
    .drop("genotypes_called")\
    .drop("called")


def checkSize(df, col2Check, colName, default = 0):
    '''
    Check if alt value for array exists and return value otherwise default value is returned
    '''
    return \
    df.withColumn(colName, f.when(f.size(f.col(col2Check)) > 1, f.col(col2Check)[1]).otherwise(default))

def extractValues(df, fields_to_include = []):
    fields_to_select = [ "contigName", "start", "end", "referenceAllele", "alternateAlleles", "filters", "qual", "INFO_OLD_MULTIALLELIC AS multiAllelic",\
    "uid", "CHR", "POS", "REF", "ALT", "AC", "aaf", "nHomozygous[0] AS homRR", "nHet AS hetRA", "homAA", "pValueHwe", "callRate", "nUncalled AS noCall",\
    "depth_quartiles", "het_alleleDepth_quartiles", "hom_alleleDepth_quartiles", "genotypes" ]
    fields_to_select.extend(fields_to_include)
    
    '''
    Run function as: pvcf_df.transform(lambda df: toPvcfFormat(df, ["hg19_uid AS INFO_hg19_uid", "INFO_OLD_MULTIALLELIC AS INFO_oldMultiAllelic"]))
    This function will calculate min, med and max het, hom AB ratios, depths. Convert array alleleCounts, allelefrequencies and nHomozygous counts to aaf and homAA
    '''
    
    return \
    df\
    .transform(lambda df: checkSize(df, "nHomozygous", "homAA"))\
    .transform(lambda df: checkSize(df, "alleleCounts", "AC"))\
    .transform(lambda df: checkSize(df, "alleleFrequencies", "aaf"))\
    .withColumn("depths", f.expr("transform(genotypes, g -> if(!array_contains(g.calls, -1), g.depth, -1))"))\
    .withColumn("depths", f.array_sort(f.expr("filter(depths, d -> d > -1)")))\
    .withColumn("depths_n", f.size(f.col("depths")))\
    .withColumn("depth_quartiles", f.array(f.col("depths")[0], f.col("depths")[(f.col("depths_n") *  0.5).cast("Int")], f.col("depths")[f.col("depths_n") - 1]))\
    .withColumn("het_alleleDepths", \
    f.expr("transform(genotypes, g -> \
    if(array_contains(g.calls, 0) AND array_contains(g.calls, 1), g.alleleDepths[1] / (g.alleleDepths[0] + g.alleleDepths[1]), -1))"))\
    .withColumn("het_alleleDepths", f.array_sort(f.expr("filter(het_alleleDepths, a -> a > -1)")))\
    .withColumn("het_alleleDepth_quartiles", f.when(f.col("nHet") > 0,\
    f.struct(f.col("het_alleleDepths")[0].alias("min"), f.col("het_alleleDepths")[(f.col("nHet") * 0.5).cast("Int")].alias("median"), \
    f.col("het_alleleDepths")[f.col("nHet") - 1].alias("max")))\
    .otherwise(f.struct(f.lit(0).alias("min"),f.lit(0).alias("median"),f.lit(0).alias("max"))))\
    .withColumn("hom_alleleDepths", \
    f.expr("transform(genotypes, g -> if(g.calls[0] == 1 AND g.calls[1] == 1, g.alleleDepths[1] / (g.alleleDepths[0] + g.alleleDepths[1]), -1))"))\
    .withColumn("hom_alleleDepths", f.array_sort(f.expr("filter(hom_alleleDepths, a -> a > -1)")))\
    .withColumn("hom_alleleDepth_quartiles", f.when(f.col("homAA") > 0,\
    f.struct(f.col("hom_alleleDepths")[0].alias("min"), f.col("hom_alleleDepths")[(f.col("homAA") * 0.5).cast("Int")].alias("median"),\
    f.col("hom_alleleDepths")[f.col("homAA")].alias("max")))\
    .otherwise(f.struct(f.lit(0).alias("min"),f.lit(0).alias("median"),f.lit(0).alias("max"))))\
    .selectExpr(fields_to_select)

