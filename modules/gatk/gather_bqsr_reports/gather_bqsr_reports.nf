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


process GATK_GATHER_BQSR_REPORTS {
    tag "${sampleId}"

    container "quay.io/seqeralabs/gatk4-germline-snps-indels"
    memory "16 GB"
    cpus 16

    input:
    tuple val(sampleId), path(input_bqsr_reports)

    output:
    tuple val(sampleId),
            path("${sampleId}.recal_data.csv")

    script:
    input_bqsr_params = input_bqsr_reports.join(" -I ")

    """
    ${params.gatk_path} --java-options ${params.java_opts} \
                        GatherBQSRReports \
                        -I ${input_bqsr_params} \
                        -O "${sampleId}.recal_data.csv"
    """

    stub:

    """
    touch "${sampleId}.recal_data.csv"
    """
}
