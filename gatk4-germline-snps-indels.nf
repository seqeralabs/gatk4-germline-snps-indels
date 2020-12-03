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
// Derive file names and location from the params.yaml
//================================================================================

fastq_files_list = file(params.input_fofn)

//================================================================================
// Include sub-workflows and (soft) override workflow-level parameters
//================================================================================

include { FORMAT_CONVERSION } from "./workflows/format_conversion/format_conversion.nf"
include { PREPROCESSING_MAPPING } from "./workflows/preprocessing_mapping/preprocessing_mapping.nf"
include { QUALITY_RECALIBRATION } from "./workflows/quality_recalibration/quality_recalibration.nf"
include { VARIANT_DISCOVERY } from "./workflows/variant_discovery/variant_discovery.nf"


//================================================================================
// Prepare channels
//================================================================================


// Convert input manifest to a channel.
fastq_params_ch = channel.fromPath(fastq_files_list)
        .splitText(keepHeader: false)
        .map { line ->
            cols = line.tokenize('\t')
            [
                    file(cols[2]), // fastq_1
                    file(cols[3]), // fastq_2
                    cols[6], // run_date
                    cols[1], // sample_name
                    cols[4], // library_name
                    cols[7], // platform_name
                    cols[5], // platform_name
                    cols[0], // readgroup_name
                    cols[8] // sequencing_center

            ]
        }

//================================================================================
// Main workflow
//================================================================================

workflow {

    FORMAT_CONVERSION({ fastq_params_ch })

    PREPROCESSING_MAPPING(FORMAT_CONVERSION.out)

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

