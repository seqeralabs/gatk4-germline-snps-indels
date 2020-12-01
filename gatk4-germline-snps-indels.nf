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


//================================================================================
// Read and derive file names and location from the params
//================================================================================

unmapped_bams = file(params.unmapped_bams_list)

//================================================================================
// Include modules and (soft) override sub-workflow-level parameters
//================================================================================

include { FORMAT_CONVERSION } from "./workflows/format_conversion/format_conversion.nf"
include { PREPROCESSING_MAPPING } from "./workflows/preprocessing_mapping/preprocessing_mapping.nf"
include { QUALITY_RECALIBRATION } from "./workflows/quality_recalibration/quality_recalibration.nf"
include { VARIANT_DISCOVERY } from "./workflows/variant_discovery/variant_discovery.nf"


//================================================================================
// Prepare channels
//================================================================================

/*

If we use a tab-delimited files with sampleId in one column and path to unaligned bam in the second,
it will eliminate the need make assumptions about base filename structure

*/
unmapped_bams_ch = channel.fromPath(unmapped_bams)
        .splitText()
        .map { line -> [line.tokenize("\t")[0], file(line.tokenize("\t")[1].trim())] }


//================================================================================
// Main workflow
//================================================================================

workflow {

    PREPROCESSING_MAPPING(unmapped_bams_ch)

    QUALITY_RECALIBRATION(PREPROCESSING_MAPPING.out)

    VARIANT_DISCOVERY(QUALITY_RECALIBRATION.out)

    VARIANT_DISCOVERY.out
            .map {
                sampleId, vcfFile -> "${sampleId}\ts3://${vcfFile}"
            }
            .collectFile(
                    name: 'merged_vcfs.tsv', newLine: true, storeDir: "${params.outdir}"
            )
}

