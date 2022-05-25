vcf_name=${1%.vcf*}
echo $vcf_name

java -Xmx32g -jar $picard LiftoverVcf I=$1 O=${vcf_name}.hg38.vcf.gz R=/ref_genomes/hg38.nochr.fa CHAIN=/hg19ToHg38.over.chain_nochr.gz REJECT=${vcf_name}.failed.vcf WRITE_ORIGINAL_POSITION=true RECOVER_SWAPPED_REF_ALT=true WRITE_ORIGINAL_ALLELES=true ALLOW_MISSING_FIELDS_IN_HEADER=true TMP_DIR=. MAX_RECORDS_IN_RAM=2500000
