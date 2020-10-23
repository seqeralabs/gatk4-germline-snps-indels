nextflow.enable.dsl = 2

process GATK_COLLECT_VARIANT_CALLING_METRICS {

    input:
    path(input_vcf)
    path(input_vcf_index)
    path(dbsnp_vcf)
    path(dbsnp_vcf_index)
    path(interval_list)
    path(ref_dict)

    output:
    path("${params.metrics_file_prefix}.variant_calling_detail_metrics")
    path("${params.metrics_file_prefix}.variant_calling_summary_metrics")

    script:
    """
    set -euo pipefail

    gatk --java-options -Xms6g CollectVariantCallingMetrics \
                               --INPUT ${input_vcf} \
                               --DBSNP ${dbsnp_vcf} \
                               --SEQUENCE_DICTIONARY ${ref_dict} \
                               --OUTPUT ${params.metrics_filename_prefix} \
                               --THREAD_COUNT 8 \
                               --TARGET_INTERVALS ${interval_list}

    """
}
