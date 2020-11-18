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

process GATK_MARK_DUPLICATES {
    tag "${sampleId}"

    container "quay.io/seqeralabs/gatk4-germline-snps-indels"
    memory "32 GB"
    cpus 16

    input:
    val(sampleId)
    path(input_mapped_merged_bam)

    output:
    val(sampleId)
    path("*_merged.deduped.bam")
    path("*_merged.deduped.metrics.txt")

    script:

    """
    ${params.gatk_path} --java-options "-Dsamjdk.compression_level=${params.compression_level} ${params.java_opts}" \
                        MarkDuplicates \
                        --INPUT ${input_mapped_merged_bam} \
                        --OUTPUT ${sampleId}_merged.deduped.bam \
                        --METRICS_FILE ${sampleId}_merged.deduped.metrics.txt \
                        --VALIDATION_STRINGENCY SILENT \
                        --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
                        --ASSUME_SORT_ORDER "queryname" \
                        --CREATE_MD5_FILE true
    """

}
