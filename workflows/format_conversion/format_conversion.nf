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
// Include modules and (soft) override module-level parameters
//================================================================================

include { GATK_PAIRED_FASTQ_TO_UNMAPPED_BAM } from "../../modules/gatk/paired_fastq_to_unmapped_bam/paired_fastq_to_unmapped_bam.nf"


//================================================================================
// Main workflow
//================================================================================

workflow FORMAT_CONVERSION {
    take:
    fastq_files

    main:
    GATK_PAIRED_FASTQ_TO_UNMAPPED_BAM(fastq_files)

    emit:
    GATK_PAIRED_FASTQ_TO_UNMAPPED_BAM.out
}

