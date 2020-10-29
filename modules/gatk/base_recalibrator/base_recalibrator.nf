nextflow.enable.dsl = 2

params.container = "broadinstitute/gatk:4.1.8.1"
params.gatk_path = "/gatk/gatk"
params.memory = '16'
params.cpus = 16
params.java_opts = ""


process GATK_BASE_RECALIBRATOR {
    tag "${sampleId}_${subgroup_name}"

    container params.container
    memory "${params.memory} GB"
    cpus params.cpus

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
    subgroup_trimmed = subgroup.trim().split("\t").join(" -L ")

    """
    ${params.gatk_path} --java-options ${params.java_opts} \
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

workflow test {

}
