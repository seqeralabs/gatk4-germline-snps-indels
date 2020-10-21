# Download the data using gsutil and the official gatk-test-data buckets on Google cloud.

# - Install gsutil https://cloud.google.com/storage/docs/gsutil_install#linux
# - Login to your profile
# - Run this script
# - (Profit!)

set -uex

# INPUT_BAM
gsutil cp gs://gatk-test-data/wgs_bam/NA12878_24RG_hg38/NA12878_24RG_small.hg38.bam ./
gsutil cp gs://gatk-test-data/wgs_bam/NA12878_24RG_hg38/NA12878_24RG_small.hg38.bai ./

# REFERENCE FILES
gsutil cp gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta ./
gsutil cp gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai ./
gsutil cp gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dict ./

# INTERVALS
gsutil cp gs://gatk-test-data/intervals/hg38_wgs_scattered_calling_intervals.txt ./
