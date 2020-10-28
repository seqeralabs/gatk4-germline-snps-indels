nextflow.enable.dsl = 2

params.container = "broadinstitute/gatk:4.1.8.1"
params.gatk_path = "/gatk/gatk"
params.memory = '16'
params.cpus = 16
params.java_opts = "-XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10"
params.contamination = 0

process GATK_HAPLOTYPE_CALLER {
    tag "${sampleId}_${interval_chunk_name}"

    container params.container
    memory "${params.memory} GB"
    cpus params.cpus


    input:

    tuple val(sampleId),
            path(input_recal_merged_bam),
            path(input_recal_merged_bai),
            path(input_recal_merged_md5),
            val(scatter_id),
            val(interval_chunk_name),
            path(interval_list_file)

    path(ref_dict)
    path(ref_fasta)
    path(ref_fasta_fai)


    output:

    tuple val(sampleId),
            path("${sampleId}.${scatter_id.toString().padLeft(2, '0')}.${interval_chunk_name}.vcf"),
            path("${sampleId}.${scatter_id.toString().padLeft(2, '0')}.${interval_chunk_name}.vcf.idx")


    script:

    """
    set -e

    ${params.gatk_path} --java-options " -Xmx${params.memory}G ${params.java_opts}" \
                        HaplotypeCaller \
                        -R ${ref_fasta} \
                        -I ${input_recal_merged_bam} \
                        --output "${sampleId}.${scatter_id.toString().padLeft(2, '0')}.${interval_chunk_name}.vcf" \
                        -contamination ${params.contamination} \
                        -ERC GVCF \
                        -L ${interval_list_file}
    """
}

