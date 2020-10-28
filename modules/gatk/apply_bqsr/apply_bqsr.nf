
process applyBQSR {
    tag "${sampleId}_${subgroup_unmapped_name}"

    memory '16 GB'
    cpus 16

    container "broadinstitute/gatk:4.1.8.1"

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

    gatk_path = "/gatk/gatk"
    subgroup_unmapped_trimmed = subgroup_unmapped.trim().split("\t").join(" -L ")

    """
    ${gatk_path} --java-options "${params.java_opt_applybqsr}" \
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
