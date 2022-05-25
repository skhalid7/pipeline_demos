# Harmonize GWAS summstats

Given an input file (input_meta.csv) with parameters, will QC and harmonize GWAS summary stats

QC checks:
1. Filter out rare < 0.1% aaf variants
2. Filter out MAC < 25 carrier cases
3. InfoScore > 0.3

Harmonizations:
1. Set ref and alt alleles to GRCh38 reference
2. lift over to GRCh38
3. Remove I/D alleles
4. Left align indels using bcftools
