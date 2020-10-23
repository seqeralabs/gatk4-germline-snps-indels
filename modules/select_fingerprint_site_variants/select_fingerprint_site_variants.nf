nextflow.enable.dsl = 2

params.base_output_name

process GATK_SELECT_FINGERPRINT_SITE_VARIANTS {
    container "us.gcr.io/broad-gatk/gatk:4.1.1.0"

    input:
    path(input_vcf)
    path(haplotype_database)

    output:
    path("${params.base_output_name}.vcf.gz")
    path("${params.base_output_name}.vcf.gz.tbi")

    shell:

    '''
    set -euo pipefail

    function hdb_to_interval_list() {
        input=$1
        awk 'BEGIN{IFS="\t";OFS="\t";} $0~"^@"{print;next;} $0~"#CHROM"{next;} {print $1,$2,$2,"+","interval-"NR}' $1
    }

    hdb_to_interval_list !{haplotype_database} > hdb.interval_list

    gatk --java-options -Xms6g SelectVariants \
                        --variant !{input_vcf} \
                        --intervals hdb.interval_list \
                        --output !{params.base_output_name}.vcf.gz
    '''
}
