# Download the data using gsutil and the official gatk-test-data buckets on Google cloud.

# - Install gsutil https://cloud.google.com/storage/docs/gsutil_install#linux
# - Login to your profile
# - Run this script
# - (Profit!)

set -uex

# INPUT_BAM
gsutil cp gs://gatk-test-data/wgs_bam/NA12878_24RG_hg38/NA12878_24RG_small.hg38.bam ./
gsutil cp gs://gatk-test-data/wgs_bam/NA12878_24RG_hg38/NA12878_24RG_small.hg38.bai ./

# INPUT_CRAM
gsutil cp gs://gatk-test-data/wgs_cram/NA12878_20k_hg38/NA12878.cram ./
gsutil cp gs://gatk-test-data/wgs_cram/NA12878_20k_hg38/NA12878.crai ./

# REFERENCE FILES
gsutil cp gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta ./
gsutil cp gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai ./
gsutil cp gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dict ./

# INTERVALS
gsutil cp gs://gatk-test-data/intervals/test-intervals.hg38.list ./
gsutil cp gs://gatk-test-data/intervals/hg38_wgs_scattered_calling_intervals.txt ./
gsutil cp gs://gcp-public-data--broad-references/hg38/v0/scattered_calling_intervals/ ./
