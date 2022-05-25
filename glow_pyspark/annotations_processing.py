import json
import pyspark.sql.functions as f
import pandas as pd
from pyspark.sql.types import StringType, IntegerType, DoubleType, StructField, BooleanType, ArrayType, StructType, StructField

'''
Functions used to Process VEP output from GLOW 0.6. 
These functions will be used to process a dataframe with variant information into an annotated file for all protein coding Transcripts
Notebooks are run asynchronously and in the background. Notebook must have an interpreter binded to it.
'''

vep_command = [ "/ensembl-vep-release-101/vep", "--offline", "--cache", "--merged", "--format", "vcf", "--dir_cache", "/", "--output_file", "STDOUT",\
"--no_stats", "--mane", "--canonical", "--hgvs", "--shift_hgvs", "0", "--fasta", "/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz", "--vcf" ]

csq_schema = "array<struct<consequence:string,EnsemblGeneId:string,GeneName:string,exon:string,intron:string,transcript:string,hgvsc:string,\
hgvsp:string,canonical:string,mane:string,cds_position:string,protein_position:string>>"

effects_schema = "array<struct<priority:int,variantEffect:string,islof:boolean>>"

all_features = spark.read.parquet("all_features_start_end_parquet")
'''
    all_features file contains start and stops of all features in Ensembl GTF file. Format:
    CHR featureStart featureEnd geneId GeneName biotype
'''

dbnsfp = spark.read.parquet("dbNSFP_annotations")\
.withColumn("uid", f.concat_ws(":", f.col("CHR"), f.col("POS"), f.col("REF"), f.col("ALT")))\
.withColumnRenamed("rsid_151", "rsid").withColumnRenamed("", "delMissenseCount")
'''
    Format of dbnsfp file
    CHR POS REF ALT EnsemblGeneId ancestralAllele rsid missensePredictions delMissenseCount MetaSVM Revel 
    missensePredictions: predictions from 6 missense algorithms FATHMM, SIFT, Polyphen2 HDIV & HVAR, LRT, MutationAssesor 
    deleteriousMissenseCount: How many of the 6 algorithms predict the variant to be deleterious
'''

Ts2Keep = spark.read.parquet("ensembl101_transcripts")\
.withColumnRenamed("start", "cdsStart").withColumnRenamed("end", "cdsEnd") 
'''
    Ts2Keep file contains all protein coding transcripts with annotated start and stop codons
    Transcript GeneName EnsemblGeneId Strand CHR Start End CDSLength
'''

def parseCSQ(CSQ):
    '''
    This udf will parse the raw csq string to a struct
    '''
    struct = "'consequence': '{}', 'EnsemblGeneId': '{}', 'GeneName': '{}', 'exon': '{}', 'intron': '{}', 'transcript': '{}', 'hgvsc': '{}', 'hgvsp': '{}', 'canonical': '{}', 'mane': '{}',\
    'cds_position': '{}', 'protein_position': '{}'"
    csqs = CSQ.split(",")
    return "[ " + ",".join([ "{" + struct.format(csq.split("|")[1], csq.split("|")[4], csq.split("|")[3], csq.split("|")[8], csq.split("|")[9], csq.split("|")[6],\
    csq.split("|")[10], csq.split("|")[11], csq.split("|")[23], csq.split("|")[24], csq.split("|")[13], csq.split("|")[14]) + "}" for csq in csqs ]) + " ]"

parseCSQ_udf = f.udf(parseCSQ, StringType())

def harmonizeConsequence (effects: pd.Series) -> pd.Series:
    '''
    harmonize VEP effects to standardized dictionary. If > 1 effect take the most deleterious one
    effect name -> (effect priority, harmonized effect, islof)
    note: islof is deprecated
    '''
    
    effects_to_keep = { "stop_gained" : "{'priority': 1, 'variantEffect': 'stop_gained', 'islof': true}",\
    "frameshift_variant" : "{'priority': 2, 'variantEffect': 'frameshift', 'islof': true}",\
    "splice_donor_variant" : "{'priority': 3, 'variantEffect': 'splice_donor', 'islof': true}",\
    "splice_acceptor_variant" : "{'priority': 4, 'variantEffect': 'splice_acceptor', 'islof': true}",\
    "stop_lost" : "{'priority': 5, 'variantEffect': 'stop_lost', 'islof': true}",\
    "start_lost" : "{'priority': 6, 'variantEffect': 'start_lost', 'islof': true}",\
    "missense_variant" : "{'priority': 7, 'variantEffect': 'missense', 'islof': false}",\
    "inframe_insertion" : "{'priority': 8, 'variantEffect': 'inframe_indel', 'islof': false}",\
    "inframe_deletion" : "{'priority': 8, 'variantEffect': 'inframe_indel', 'islof': false}",\
    "splice_region_variant" : "{'priority': 9, 'variantEffect': 'splice_region', 'islof': false}",\
    "synonymous_variant" : "{'priority': 10, 'variantEffect': 'synonymous', 'islof': false}",\
    "start_retained_variant" : "{'priority': 10, 'variantEffect': 'synonymous', 'islof': false}",\
    "stop_retained_variant" : "{'priority': 10, 'variantEffect': 'synonymous', 'islof': false}",\
    "5_prime_UTR_variant" : "{'priority': 11, 'variantEffect': '5_prime_UTR', 'islof': false}",\
    "3_prime_UTR_variant" : "{'priority': 12, 'variantEffect': '3_prime_UTR', 'islof': false}",\
    "intron_variant" : "{'priority': 13, 'variantEffect': 'intronic', 'islof': false}",\
    "upstream_gene_variant" : "{'priority': 14, 'variantEffect': 'upstream', 'islof': false}",\
    "downstream_gene_variant" : "{'priority': 15, 'variantEffect': 'downstream', 'islof': false}" }

    return effects.apply(lambda consequences: "[ " + ",".join([ effects_to_keep[effect] if effect in effects_to_keep else "{'priority': 16, 'variantEffect': 'other', 'islof': false}" for effect in consequences]) + " ]")

harmonizeConsequence_udf = f.pandas_udf(harmonizeConsequence, returnType=StringType())

def getNearestGenes(df, all_features_df, cols_2_include = []):
    '''
    Join with all_features
    Group By variant and get all features within +/- 500kbps
    For each variant identify 4 nearest protein coding genes and 4 nearest non-protein coding features
    After groupBy original columns get lost, to keep those add them in cols2_include
    '''
    cols_2_groupBy = [ "uid", "CHR", "ref", "alt", "pos" ]
    cols_2_groupBy.extend(cols_2_include)
    
    return df.join(f.broadcast(all_features_df), ["CHR"])\
    .filter((f.col("pos").between(f.col("featureStart"), f.col("featureEnd"))) | \
    (f.col("pos").between(f.col("featureEnd"), f.col("featureEnd") + 500000)) | \
    (f.col("pos").between(f.col("featureStart") - 500000, f.col("featureStart"))))\
    .withColumn("distance", f.when((f.col("pos").between(f.col("featureStart"), f.col("featureEnd"))), f.lit(0).cast("Int"))\
    .when((f.col("pos").between(f.col("featureEnd"), f.col("featureEnd") + 500000)), f.col("featureEnd") - f.col("pos"))\
    .when((f.col("pos").between(f.col("featureStart") - 500000, f.col("featureStart"))), f.col("featureStart") - f.col("pos")))\
    .withColumn("biotypeCategory", f.when(f.col("biotype") == "protein_coding", "protein_coding").otherwise("other"))\
    .withColumn("absDistance", f.abs(f.col("distance")))\
    .groupBy(cols_2_groupBy).pivot("biotypeCategory", ["protein_coding", "other"])\
    .agg(f.array_sort(f.collect_set(f.struct("absDistance", "GeneName", f.col("geneId").alias("EnsemblGeneId"), "biotype", "distance"))).alias("nearestGenes"))\
    .withColumn("nearestProteinCodingGenes", f.slice("protein_coding", 1, 4))\
    .withColumn("nearestNonCodingFeatures", f.slice("other", 1, 4))\
    .drop("protein_coding").drop("other")

def createLofFlags(df):
    '''
    This funciton will create the lof Flags: 5_UTR_splice, 3_UTR_splice, GC_2_GT_donor and isEndTrunc, ancestral_allele
    Identify splice donor/acceptor variants in 3' and 5' UTRs, splice donor variants mutations GC -> GTs, variants in the last 5%
    lofFilters: comma seperated string of all flags. "" if no flag.
    '''

    return df\
    .withColumn("ancestralLof", f.when((f.col("isAncestralAllele")) & (f.col("effects.priority") >= 4), "ancestral_allele").otherwise("NA"))\
    .withColumn("isUTRSplice", f.when(((f.col("effects.variantEffect") == "splice_donor") | \
    (f.col("effects.variantEffect") == "splice_acceptor")) & ((f.col("POS") <= f.col("cdsStart")) & (f.col("strand") == "+")), "5_UTR_splice")\
    .when(((f.col("effects.variantEffect") == "splice_donor") | (f.col("effects.variantEffect") == "splice_acceptor")) & \
    ((f.col("POS") >= f.col("cdsEnd")) & (f.col("strand") == "+")), "3_UTR_splice")\
    .when(((f.col("effects.variantEffect") == "splice_donor") | (f.col("effects.variantEffect") == "splice_acceptor")) & \
    ((f.col("POS") >= f.col("cdsEnd")) & (f.col("strand") == "-")), "5_UTR_splice")\
    .when(((f.col("effects.variantEffect") == "splice_donor") | (f.col("effects.variantEffect") == "splice_acceptor")) & \
    ((f.col("POS") <= f.col("cdsStart")) & (f.col("strand") == "-")), "3_UTR_splice")\
    .otherwise("NA"))\
    .withColumn("isGC2GT", f.when((f.col("effects.variantEffect") == "splice_donor") & (f.col("HGVSc").contains("+2C>T")), "GC_2_GT_donor").otherwise("NA"))\
    .withColumn("isEndTrunc", f.when(((f.col("effects.variantEffect") == "frameshift") | (f.col("effects.variantEffect") == "stop_gained")) & \
    (f.col("CDS_Position") / f.col("cdsLength") >= 0.95), "EndTrunc").otherwise("NA"))\
    .withColumn("lofFilters", f.array_remove(f.array("isUTRSplice", "isGC2GT", "isEndTrunc", "ancestralLof"), "NA"))\
    .withColumn("lofFilters", f.when(f.size(f.col("lofFilters")) > 0, f.concat_ws(",", f.col("lofFilters"))).otherwise("NA"))

def mergeEffects(df, cols_2_include = []):

    '''
    Merge all effects for a variant per gene -->  gene - variant
    Create a column for variant effect on canonical transcript
    Create a column for most deleterious consequence
    '''

    cols_2_groupBy = [ "uid", "CHR", "ref", "alt", "pos", "rsid", "isAncestralAllele", "GeneName", "EnsemblGeneId" ]
    cols_2_groupBy.extend(cols_2_include)

    return df.withColumn("affectsCanonical", f.expr("filter(completeTranscriptAnnotations, c -> c.canonical = 'yes')")[0].isNotNull())\
    .groupBy(cols_2_groupBy)\
    .agg(f.array_sort(f.collect_list(f.struct("priority", "variantEffect", "delMissenseCount", "missensePredictions", "revel", "completeTranscriptAnnotations"))).alias("allEffects"),\
    f.first(f.when(f.col("affectsCanonical"), f.struct("priority", "variantEffect", "delMissenseCount", "missensePredictions", "revel", "completeTranscriptAnnotations"))).alias("canonicalAnnotations"))\
    .withColumn("canonicalAnnotations", f.struct("canonicalAnnotations.priority", "canonicalAnnotations.variantEffect", "canonicalAnnotations.delMissenseCount", "canonicalAnnotations.missensePredictions",\
    "canonicalAnnotations.revel", f.expr("filter(canonicalAnnotations.completeTranscriptAnnotations, c -> c.canonical == 'yes')")[0].alias("completeTranscriptAnnotations")))\
    .withColumn("mostDeleteriousEffect", f.col("allEffects")[0])

def getMostDeleteriousAnnotation(df):

    '''
    Parse mostDeleteriousEffect column
    Definition of HC lof: None of the transcripts are flagged by lofFlags
    '''

    return df\
    .withColumn("lofFlags", f.size(f.expr("filter(mostDeleteriousEffect.completeTranscriptAnnotations.lofFilters, lf -> lf != 'NA')")) > 0)\
    .withColumn("HClof", f.when((f.col("lofFlags")) & (f.col("mostDeleteriousEffect.priority") <= 4), "LC")\
    .when(~(f.col("lofFlags")) & (f.col("mostDeleteriousEffect.priority") <= 4), "HC")\
    .when((f.col("mostDeleteriousEffect.priority") == 5) | (f.col("mostDeleteriousEffect.priority") == 6), "LC")\
    .otherwise("NA"))\
    .withColumn("mostDeleteriousEffect", f.struct("mostDeleteriousEffect.variantEffect", f.col("mostDeleteriousEffect.priority").cast("String").alias("priority"),\
    f.concat_ws(":", "mostDeleteriousEffect.completeTranscriptAnnotations.Transcript").alias("transcripts"), f.array_contains("mostDeleteriousEffect.completeTranscriptAnnotations.canonical", "yes").alias("affectsCanonical"),\
    f.col("mostDeleteriousEffect.delMissenseCount").cast("String").alias("delMissenseCount"), "mostDeleteriousEffect.missensePredictions", f.col("mostDeleteriousEffect.revel").cast("String").alias("revel"),\
    f.concat_ws(":", "mostDeleteriousEffect.completeTranscriptAnnotations.hgvsc").alias("hgvsc"),\
    f.concat_ws(":", "mostDeleteriousEffect.completeTranscriptAnnotations.hgvsp").alias("hgvsp"), f.concat_ws(":", "mostDeleteriousEffect.completeTranscriptAnnotations.cds").alias("cds"),\
    f.concat_ws(":", "mostDeleteriousEffect.completeTranscriptAnnotations.aa").alias("aa"), f.concat_ws(":", "mostDeleteriousEffect.completeTranscriptAnnotations.exon").alias("exon"),\
    f.concat_ws(":", "mostDeleteriousEffect.completeTranscriptAnnotations.intron").alias("intron"),\
    f.concat_ws(":", "mostDeleteriousEffect.completeTranscriptAnnotations.lofFilters").alias("lofFilters"), "HClof"))\
    .drop("lofFlags").drop("HClof")

def getCanonicalAnnotations(df):

    '''
    Get Canonical Annotation and identify if it's HC or not.
    '''

    return df\
    .withColumn("HClof", f.when(f.col("canonicalAnnotations.variantEffect").isNull(), None).when((f.col("canonicalAnnotations.completeTranscriptAnnotations.lofFilters") != "NA") & (f.col("canonicalAnnotations.priority") <= 4), "LC")\
    .when((f.col("canonicalAnnotations.completeTranscriptAnnotations.lofFilters") == "NA") & (f.col("canonicalAnnotations.priority") <= 4), "HC")\
    .when((f.col("canonicalAnnotations.priority") == 5) | (f.col("canonicalAnnotations.priority") == 6), "LC")\
    .otherwise("NA"))\
    .withColumn("lofFilters", f.when(f.col("canonicalAnnotations.variantEffect").isNull(), None).when(f.col("canonicalAnnotations.completeTranscriptAnnotations.lofFilters") == "NA", "NA")\
    .otherwise(f.concat_ws(",", f.col("canonicalAnnotations.completeTranscriptAnnotations.lofFilters"))))\
    .withColumn("canonicalEffect", f.when(f.col("canonicalAnnotations.variantEffect").isNotNull(), f.struct("canonicalAnnotations.variantEffect", f.col("canonicalAnnotations.priority").cast("String"),\
    f.col("canonicalAnnotations.completeTranscriptAnnotations.Transcript").alias("transcripts"),\
    f.lit(True).cast("String").alias("affectsCanonical"), f.col("canonicalAnnotations.delMissenseCount").cast("String").alias("delMissenseCount"), "canonicalAnnotations.missensePredictions",\
    f.col("canonicalAnnotations.revel").cast("String").alias("revel"),\
    "canonicalAnnotations.completeTranscriptAnnotations.hgvsc", "canonicalAnnotations.completeTranscriptAnnotations.hgvsp", "canonicalAnnotations.completeTranscriptAnnotations.cds",\
    "canonicalAnnotations.completeTranscriptAnnotations.aa", "canonicalAnnotations.completeTranscriptAnnotations.exon", "canonicalAnnotations.completeTranscriptAnnotations.intron", "lofFilters", "HClof")).otherwise(None))\
    .drop("HClof").drop("canonicalAnnotations").drop("lofFilters")

def getCompleteAnnotations(df, csq_col, dbnsfp_df, Ts2Keep_df, cols_2_include = []):
    '''
    Parse CSQ field (vcf must be read in with option("flattenInfoFields", False) and select relevant fields
    Merge with dbnsfp annotations containing Missense predictions
    Join with all protein coding transcripts with annotated start and stop codons
    Harmonize consequences and assign effect priority
    Null out empty strings

    Generate lof flags with createLofFlags function
    For each variant effect merge all transcripts with the same effect

    Definition of HC -> stop_gained, frameshift, splice_donor, splice_acceptor without any lofFlags
    Definition of LC -> stop_gained, frameshift, splice_donor, splice_acceptor with lofFlags and all start lost, stop lost variants

    After groupBy original columns get lost, to keep those add them in cols2_include

    Returns variant idenitifiers, GeneName, EnsemblGeneId, canonicalEffects, mostDeleteriousEffect, allEffects (an array containing effect on all transcripts)

    '''

    dbnsfp_annotations = [ "isAncestralAllele", "delMissenseCount", "missensePredictions", "revel" ]
    cols_2_groupBy = [ "uid", "CHR", "ref", "alt", "pos", "rsid", "GeneName", "EnsemblGeneId", "effects.variantEffect", "effects.priority" ] + dbnsfp_annotations
    cols_2_groupBy.extend(cols_2_include)


    return df\
    .withColumn("INFO_CSQ", f.from_json(parseCSQ_udf(f.col(csq_col)), csq_schema))\
    .withColumn("INFO_CSQ", f.explode(f.col("INFO_CSQ")))\
    .select("*", "INFO_CSQ.*")\
    .drop("INFO_CSQ").drop("attributes").drop("val").drop(csq_col)\
    .join(dbnsfp_df, ["uid", "CHR", "POS", "REF", "ALT", "EnsemblGeneId"], "left_outer")\
    .join(f.broadcast(Ts2Keep_df), ["EnsemblGeneId", "GeneName", "Transcript", "CHR"])\
    .withColumn("effects", f.from_json(harmonizeConsequence_udf(f.split(f.col("consequence"), "&")), effects_schema))\
    .withColumn("effects", f.array_sort("effects")[0])\
    .withColumn("consequence", f.concat_ws(",", f.col("consequence")))\
    .withColumn("CDS_position", f.split(f.col("CDS_position"), "-")[0])\
    .withColumn("Protein_position", f.split(f.col("Protein_position"), "-")[0])\
    .withColumn("cds", f.when(f.length(f.col("CDS_Position")) > 0, f.concat_ws("_", f.col("CDS_Position"), f.col("cdsLength"))).otherwise(None))\
    .withColumn("aa", f.when(f.length(f.col("Protein_position")) > 0, f.concat_ws("_", f.col("Protein_position"), (f.col("cdsLength") / 3).cast("Int"))).otherwise(None))\
    .withColumn("exon", f.when(f.length(f.col("exon")) > 0, f.split(f.col("exon"), "/")).otherwise(None))\
    .withColumn("exon", f.when(f.col("exon").isNotNull(), f.concat_ws("_", f.col("exon")[0], f.col("exon")[1])).otherwise(None))\
    .withColumn("intron", f.when(f.length(f.col("intron")) > 0, f.split(f.col("intron"), "/")).otherwise(None))\
    .withColumn("intron", f.when(f.col("INTRON").isNotNull(), f.concat_ws("_", f.col("intron")[0], f.col("intron")[1])).otherwise(None))\
    .withColumn("hgvsc", f.split(f.col("hgvsc"), ":")[1])\
    .withColumn("hgvsp", f.split(f.col("hgvsp"), ":")[1])\
    .withColumn("isAncestralAllele", f.when(f.col("alt") == f.col("AncestralAllele"), True).otherwise(False))\
    .withColumn("delMissenseCount", f.when(f.col("effects.variantEffect") == "missense", f.col("delMissenseCount")).otherwise(None))\
    .withColumn("missensePredictions", f.when(f.col("effects.variantEffect") == "missense", f.col("missensePredictions")).otherwise(None))\
    .withColumn("revel", f.when(f.col("revel") == ".", None).otherwise(f.col("revel")))\
    .withColumn("canonical", f.when(f.col("canonical") == "YES", "yes").otherwise("no"))\
    .transform(lambda df : createLofFlags(df))\
    .groupBy(cols_2_groupBy)\
    .agg(f.collect_list(f.struct("Transcript", "canonical", "hgvsc", "hgvsp", "cds", "aa", "exon", "intron", "lofFilters")).alias("completeTranscriptAnnotations"))\
    .transform(lambda df: mergeEffects(df, cols_2_include))\
    .transform(lambda df: getCanonicalAnnotations(df))\
    .transform(lambda df: getMostDeleteriousAnnotation(df))

