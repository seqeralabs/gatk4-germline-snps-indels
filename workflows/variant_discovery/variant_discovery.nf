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

ref_fasta = file(params.fasta)
ref_fasta_fai = file("${params.fasta}.fai")
ref_dict = file(params.fasta.replace(".fasta", ".dict"))
scattered_calling_interval = file(params.scattered_calling_interval)

//================================================================================
// Include modules and (soft) override module-level parameters
//================================================================================

include { GATK_HAPLOTYPE_CALLER } from "../../modules/gatk/haplotype_caller/haplotype_caller.nf" addParams(params.GATK_HAPLOTYPE_CALLER)
include { GATK_MERGE_VCFS } from "../../modules/gatk/merge_vcfs/merge_vcfs.nf" addParams(params.GATK_MERGE_VCFS)


//================================================================================
// Prepare channels
//================================================================================

// Prepare channel for scattered calling intervals: input to haplotypecaller

calling_scatter_counter = 0

calling_intervals = channel.fromPath(scattered_calling_interval)
        .splitText()
        .map { line ->
            calling_scatter_counter = calling_scatter_counter + 1
            [calling_scatter_counter, line.split("/")[5], line.trim()]
        }


//================================================================================
// Main workflow
//================================================================================

workflow VARIANT_DISCOVERY {
    take:
    bam_and_internal

    main:
    GATK_HAPLOTYPE_CALLER(
            bam_and_internal.combine(calling_intervals),
            ref_dict,
            ref_fasta,
            ref_fasta_fai
    )

    GATK_MERGE_VCFS(
            GATK_HAPLOTYPE_CALLER.out.groupTuple()
    )

    emit:
    GATK_MERGE_VCFS.out


}
