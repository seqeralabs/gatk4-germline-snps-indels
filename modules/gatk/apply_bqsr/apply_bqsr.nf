nextflow.enable.dsl = 2

params.gatk_path = "/gatk/gatk"
params.java_opts = ""


process GATK_APPLY_BQSR {
    tag "${sampleId}_${subgroup_unmapped_name}"

    container = "broadinstitute/gatk:4.1.8.1"
    memory "16 GB"
    cpus 16


    input:
    tuple val(sampleId),
            path(input_mapped_merged_marked_sorted_bam),
            path(input_mapped_merged_marked_sorted_bai),
            path(input_mapped_merged_marked_sorted_md5),
            path(input_merged_bqsr_report),
            val(scatter_id),
            val(subgroup_unmapped_name),
            val(subgroup_unmapped)

    path(ref_dict)
    path(ref_fasta)
    path(ref_fasta_fai)


    output:
    tuple val(sampleId),
            path("${sampleId}.${scatter_id.toString().padLeft(2, '0')}.${subgroup_unmapped_name}.recal.bam")

    script:

    subgroup_unmapped_trimmed = subgroup_unmapped.trim().split("\t").join(" -L ")

    """
    ${params.gatk_path} --java-options "${params.java_opts}" \
                        ApplyBQSR \
                        -R ${ref_fasta}  \
                        -I ${input_mapped_merged_marked_sorted_bam}  \
                        -O "${sampleId}.${scatter_id.toString().padLeft(2, '0')}.${subgroup_unmapped_name}.recal.bam"  \
                        -L ${subgroup_unmapped_trimmed} \
                        -bqsr ${input_merged_bqsr_report} \
                        --static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30 \
                        --add-output-sam-program-record \
                        --create-output-bam-md5 \
                        --use-original-qualities
    """
}
