nextflow.enable.dsl = 2

params.gatk_path = "/gatk/gatk"
params.java_opts = ""

process GATK_PAIRED_FASTQ_TO_UNMAPPED_BAM {
    tag "${sampleId}"

    container "broadinstitute/gatk:4.1.8.1"
    memory "32 GB"
    cpus 16


    input:
    tuple file(fastq_1),
            file(fastq_2),
            val(run_date),
            val(sample_name),
            val(library_name),
            val(platform_name),
            val(platform_unit),
            val(readgroup_name),
            val(sequencing_center)

    output:
    tuple val(sample_name), path("${readgroup_name}.unmapped.bam")

    script:
    """
    sleep 60
    ${params.gatk_path} --java-options "${params.java_opts}" \
                        FastqToSam \
                        --FASTQ ${fastq_1} \
                        --FASTQ2 ${fastq_2} \
                        --OUTPUT ${readgroup_name}.unmapped.bam \
                        --READ_GROUP_NAME ${readgroup_name} \
                        --SAMPLE_NAME ${sample_name} \
                        --LIBRARY_NAME ${library_name} \
                        --PLATFORM_UNIT ${platform_unit} \
                        --RUN_DATE ${run_date} \
                        --PLATFORM ${platform_name} \
                        --SEQUENCING_CENTER ${sequencing_center}
    """
}
