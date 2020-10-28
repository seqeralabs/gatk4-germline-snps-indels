
process baseRecalibrator {
    tag "${sampleId}_${subgroup_name}"

    memory '16 GB'
    cpus 16

    container "broadinstitute/gatk:4.1.8.1"

    input:
    tuple val(sampleId),
            path(input_mapped_merged_marked_sorted_bam),
            path(input_mapped_merged_marked_sorted_bai),
            path(input_mapped_merged_marked_sorted_md5),
            val(subgroup_name),
            val(subgroup)

    path(ref_dict)
    path(ref_fasta)
    path(ref_fasta_fai)

    path(dbSNP_vcf)
    path(dbSNP_vcf_index)
    path(known_indels_mills)
    path(known_indels_mills_index)
    path(known_indels_dbSNP)
    path(known_indels_dbSNP_index)

    output:
    tuple val(sampleId),
            path("${sampleId}_recalibration_report_${subgroup_name}.recal_data.csv")

    script:
    gatk_path = "/gatk/gatk"
    subgroup_trimmed = subgroup.trim().split("\t").join(" -L ")

    """
    ${gatk_path} --java-options ${params.java_opt_baserecal} \
    BaseRecalibrator \
    -R ${ref_fasta} \
    -I ${input_mapped_merged_marked_sorted_bam} \
    --use-original-qualities \
    -O "${sampleId}_recalibration_report_${subgroup_name}.recal_data.csv" \
    --known-sites ${dbSNP_vcf} \
    --known-sites ${known_indels_mills} \
    --known-sites ${known_indels_dbSNP} \
    -L ${subgroup_trimmed}
    """
}
