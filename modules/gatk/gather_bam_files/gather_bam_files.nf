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
params.java_opts = ""
params.compression_level = 5

process GATK_GATHER_BAM_FILES {
    tag "${sampleId}"

    container "quay.io/seqeralabs/gatk4-germline-snps-indels"
    memory "16 GB"
    cpus 16

    input:
    tuple val(sampleId), path(input_recalibrated_bams)

    output:
    tuple val(sampleId),
            path("${sampleId}.recal.merged.bam"),
            path("${sampleId}.recal.merged.bai"),
            path("${sampleId}.recal.merged.bam.md5")

    script:

    inputs_bams_to_merge = input_recalibrated_bams
            .sort(false) { a, b -> a.getBaseName() <=> b.getBaseName() }.join(" --INPUT ")

    """
    ${params.gatk_path} --java-options "-Dsamjdk.compression_level=${params.compression_level}  ${params.java_opts} " \
                 GatherBamFiles \
                --INPUT ${inputs_bams_to_merge} \
                --OUTPUT "${sampleId}.recal.merged.bam" \
                --CREATE_INDEX true \
                --CREATE_MD5_FILE true
    """

    stub:

    """
    touch "${sampleId}.recal.merged.bam"
    """
}
