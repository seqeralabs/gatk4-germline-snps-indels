/*
 * Copyright (c) 2020, Seqera Labs.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 * This Source Code Form is "Incompatible With Secondary Licenses", as
 * defined by the Mozilla Public License, v. 2.0.
 */

nextflow.enable.dsl = 2

params.gatk_path = "gatk"
params.java_opts = "-Xms3000m"

process GATK_PAIRED_FASTQ_TO_UNMAPPED_BAM {
    tag "${sample_name}"

    container "quay.io/seqeralabs/gatk4-germline-snps-indels"
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

    stub:

    """
    touch ${readgroup_name}.unmapped.bam
    """
}

//================================================================================
// Module test
//================================================================================

workflow test {

    fastq_params_ch = Channel.of([
            file(params.test_read_1),
            file(params.test_read_2),
            params.run_date,
            params.sample_name,
            params.library_name,
            params.platform_name,
            params.platform_unit,
            params.readgroup_name,
            params.sequencing_center
    ])


    GATK_PAIRED_FASTQ_TO_UNMAPPED_BAM(fastq_params_ch)

}
