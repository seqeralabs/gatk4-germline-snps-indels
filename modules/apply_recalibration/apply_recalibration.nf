nextflow.enable.dsl = 2

params.recalibrated_vcf_filename
params.indel_filter_level
params.snp_filter_level
params.use_allele_specific_annotations

process GATK_APPLY_RECALIBRATION {
    container "us.gcr.io/broad-gatk/gatk:4.1.1.0"

    input:
    path(input_vcf)
    path(input_vcf_index)
    path(indels_recalibration)
    path(indels_recalibration_index)
    path(indels_tranches)
    path(snps_recalibration)
    path(snps_recalibration_index)
    path(snps_tranches)


    output:
    path("${params.recalibrated_vcf_filename}")
    path("${params.recalibrated_vcf_filename}.tbi")

    script:
    """
    set -euo pipefail

    gatk --java-options -Xms5g ApplyVQSR \
                               -O tmp.indel.recalibrated.vcf \
                               -V ${input_vcf} \
                               --recal-file ${indels_recalibration} \
                               ${params.use_allele_specific_annotations ? '--use-allele-specific-annotations' : ''} \
                               --tranches-file ${indels_tranches} \
                               --truth-sensitivity-filter-level ${indel_filter_level} \
                               --create-output-variant-index true \
                               -mode INDEL

    gatk --java-options -Xms5g ApplyVQSR \
                               -O ${recalibrated_vcf_filename} \
                               -V tmp.indel.recalibrated.vcf \
                               --recal-file ${snps_recalibration} \
                               ${params.use_allele_specific_annotations ? '--use-allele-specific-annotations' : ''} \
                               --tranches-file ${snps_tranches} \
                               --truth-sensitivity-filter-level ${snp_filter_level} \
                               --create-output-variant-index true \
                               -mode SNP
    """
}
