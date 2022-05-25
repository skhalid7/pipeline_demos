# cncd_pipelines
* [About](#about)
* [ Set up glow ](#setup)
* [pvcf to genotype lookup table](#2lookup)
* [running vep](#vep)
<a name="about"></a>
This folder contains helper functions to process pvcfs using the glow package
1. Split Multi-allelics and normalize indel representation
2. Generate het, hom, aaf, hwe_stats
3. Merge VEP annotations with dbNSFP and protein coding transcripts

Use init_scripts to set up GCP cluster

<a name="setup"></a>
## Setting up a google cloud dataproc with glow installed
In the google dataproc add the following properties:<br/>
`--properties spark:spark.hadoop.io.compression.codecs=io.projectglow.sql.util.BGZFCodec,spark:spark.jars.packages=io.projectglow:glow-spark-spark-version` <br/>
Get correct version of glow maven package from: https://glow.readthedocs.io/en/latest/getting-started.html#running-in-the-cloud<br/>
`--initialization-actions 'gs://shareef-dev/init_scripts/install_glow_py.sh'`<br/>
The init script should download the glow pywheel file (same as the glow maven package) and install it on the master node. Files such as ref-genome, vep cache should be downloaded on all nodes.
### Setting up glow & spark configs
```python
import glow
import cncd_pipelines.pvcf_transformers as pt
spark.register(glow)
```
To read in bgzipped files set the follwing cofingurations:
```
spark.conf.set("spark.hadoop.io.compression.codecs", "io.projectglow.sql.util.BGZFCodec")
```
To increase parallelization use the following setting:
```
While performing large joins update the number of shuffle partitions:
```python
spark.conf.set("spark.sql.shuffle.partitions", 1024) #update this based on number of cores available. Resulting df will have the same number of partitions as shuffle partitions
```
To run transform functions on Spark2.x monkey patch the tranform function:
```python
from pyspark.sql.dataframe import DataFrame
def transform(self, f):
    return f(self)

DataFrame.transform = transform
```
<a name="2lookup"></a>
## Converting a pvcf to a genotype look up table
**Step 1:** Read in pvcf, recode contigName (remove Chr string from Chr1 etc.), differentiate between missed and reference calls and extract calls, allele depths, GQ sccores, depth and sampleId from genotypes file. Read comments in code for exact details.

```python
input_pvcf = spark.read.format("vcf").option("flattenInfoFields", False).load(input_pvcfs_folder)\
.transform(lambda df: pt.processGenotypes(df))
```
**Step 2:** Split multi-allelic sites and correct the split multi allelic site using the correctMultiAllelics transformer. 'correctMultiAllelics' corrects split genotypes, and mimics Hails split_multi function
```python
split_pvcf = glow.transform("split_multiallelics", input_pvcf).transform(pt.correctMultiAllelics)
```
**Step 3:** Normalize variants, calcualte variant stats, remove homRR sites for the look up table and extract relevant values.
```python
split_norm_pvcf = glow.transform("normalize_variants", split_pvcf, reference_genome_path="/hg38.nochr.fa")\
.selectExpr("*", "expand_struct(call_summary_stats(genotypes))", "expand_struct(hardy_weinberg(genotypes))")\
.transform(lambda df: pt.glowFmt2uid(df))\
.withColumn("genotypes", f.expr("filter(genotypes, g -> !(g.calls[0] == 0 AND g.calls[1] == 0))"))\
.transform(lambda df: pt.extractValues(df))\
.write.parquet(path)
```
<a name="vep"></a>
## Running VEP on Glow
Set up a cluster with the install_vep.sh script. 
Use a n-8 or larger instance (don't use preemptible workers). Increase the spark executor memory and yarn over head memory: <br>
`--properties spark:spark.executor.memory=12g,spark:spark.yarn.executor.memoryOverhead=3500mb`
Important: yarn memory overhead should not be more than 25% of spark executor memory
On initialization provide the init_script: gs://shareef-dev/init_scripts/vep_installation_glow.sh
Use the following command to run vep:
```
glow.transform("pipe", site_only_vcf_df,\
cmd=vep_command, input_formatter="vcf", output_formatter="text", in_vcf_header="infer")\
.write.mode("overwrite").parquet(output_path)
```
Tips for running VEP on glow:
* Make sure the sites_df is ordered by contigName and start position
* Increase parallelization of input file by tweaking `spark.sql.files.maxPartitionBytes`
* For large pvcfs, only annotate sites and join with original pvcf later instead of passing the entire pvcf to vep
* If you run into memoryOverhead errors increase yarn memory overhead, reduce the number of executors per core

Once VEP annotated run the following block to parse annotations. This block will read in the annotated vep output as a text file. Parse out the CSQ field and format the resulting CSQ. Any site that is not within 5kb of a protein coding gene or a CSQ field will be filtered out.
```python
spark.read.parquet(annotated_file)\
filter(~f.col("text").contains("#"))\
.withColumn("splitString", f.split(f.col("text"), "\t"))\
.select(f.col("splitString")[0].alias("CHR"), f.col("splitString")[1].alias("POS"),\
f.col("splitString")[3].alias("Ref"), f.col("splitString")[4].alias("Alt"), f.col("splitString")[7].alias("INFO"))\
.withColumn("csqs", f.split(f.col("INFO"), "=")[1])\
.filter(f.col("csqs").isNotNull())\
.withColumn("uid", f.concat_ws(":", f.col("CHR"), f.col("POS"), f.col("Ref"), f.col("Alt")))\
.transform(lambda df: getCompleteAnnotations(df, "csqs", dbnsfp, Ts2Keep))\
.filter(f.col("allEffects").isNotNull())\
.write.mode("overwrite").parquet(processed_output)
```
Important: The function getCompleteAnnotations has to be loaded in seperately in another notebook
