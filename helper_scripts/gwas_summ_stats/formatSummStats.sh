output_base_path=$1
input_file=$2

PATH_2_FASTA="/ref_genomes/hg38.nochr.fa"

echo "##fileformat=VCFv4.2" > ${output_base_path}.vcf
echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" >> ${output_base_path}.vcf

summStats2vcf.py $input_file $output_base_path

cat vcf_formatted_stats >> ${output_base_path}.vcf

build=$(cat ${output_base_path}_params_file.tsv | cut -f18 | tail -1)

if [ $build == "hg19" ]
then
	echo "lift over"
	lift_over.sh ${output_base_path}.vcf
else
	echo "Already hg38"
	mv ${output_base_path}.vcf ${output_base_path}.hg38.vcf
	gzip -f ${output_base_path}.hg38.vcf
fi

swapped_alleles=$(zgrep SwappedAlleles ${output_base_path}.hg38.vcf.gz | wc -l)
echo -e "Swapped alleles: ${swapped_alleles}"

bcftools norm -f $PATH_2_FASTA -c s -O z ${output_base_path}.hg38.vcf.gz > ${output_base_path}.hg38.norm.vcf.gz

lines_2_skip=$(zcat ${output_base_path}.hg38.norm.vcf.gz | grep "##" | wc -l)

vcf2FormattedSummStats.py ${output_base_path}.hg38.norm.vcf.gz $lines_2_skip

echo -e "${output_base_path}.hg38.norm.summary_stats.formatted.tsv written"

