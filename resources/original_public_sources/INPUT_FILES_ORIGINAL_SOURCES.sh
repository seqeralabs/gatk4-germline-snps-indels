outdir - s3://nextflow-ci/gatk-ubams/results/

# From AWS public bucket
fasta - s3://broad-references/hg38/v0/Homo_sapiens_assembly38.fasta 

dbSNP_vcf - s3://broad-references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf
 
known_indels_mills - s3://broad-references/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
 
known_indels_dbSNP - s3://broad-references/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz

# Based on AWS public bucket
input_fofn - ./manifest-6samples.txt


# In this folder
sequence_grouping - ./sequence_grouping.txt

sequence_grouping_unmapped - ./sequence_grouping_with_unmapped.txt

scattered_calling_interval - ./hg38_wgs_scattered_calling_intervals.txt

