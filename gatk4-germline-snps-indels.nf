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
// Read and derive file names and location from the params.yaml
//================================================================================

unmapped_bams = file(params.unmapped_bams_list)
ref_fasta = file(params.fasta)
ref_alt = file("${params.fasta}.64.alt")
ref_amb = file("${params.fasta}.64.amb")
ref_ann = file("${params.fasta}.64.ann")
ref_bwt = file("${params.fasta}.64.bwt")
ref_pac = file("${params.fasta}.64.pac")
ref_sa = file("${params.fasta}.64.sa")
ref_fasta_fai = file("${params.fasta}.fai")
ref_dict = file(params.fasta.replace(".fasta", ".dict"))
dbSNP_vcf = file(params.dbSNP_vcf)
dbSNP_vcf_index = file("${params.dbSNP_vcf}.idx")
known_indels_mills = file(params.known_indels_mills)
known_indels_mills_index = file("${params.known_indels_mills}.tbi")
known_indels_dbSNP = file(params.known_indels_dbSNP)
known_indels_dbSNP_index = file("${params.known_indels_dbSNP}.tbi")
sequence_grouping = file(params.sequence_grouping)
sequence_grouping_unmapped = file(params.sequence_grouping_unmapped)
scattered_calling_interval = file(params.scattered_calling_interval)

//================================================================================
// Include modules and (soft) override module-level parameters
//================================================================================

include { BWA_GET_BWA_VERSION } from "./modules/bwa/get_bwa_version/get_bwa_version.nf"
include { GATK_APPLY_BQSR } from "./modules/gatk/apply_bqsr/apply_bqsr.nf" addParams(params.GATK_APPLY_BQSR)
include { GATK_BASE_RECALIBRATOR } from "./modules/gatk/base_recalibrator/base_recalibrator.nf" addParams(params.GATK_BASE_RECALIBRATOR)
include { GATK_GATHER_BAM_FILES } from "./modules/gatk/gather_bam_files/gather_bam_files.nf" addParams(params.GATK_GATHER_BAM_FILES)
include { GATK_GATHER_BQSR_REPORTS } from "./modules/gatk/gather_bqsr_reports/gather_bqsr_reports.nf" addParams(params.GATK_GATHER_BQSR_REPORTS)
include { GATK_HAPLOTYPE_CALLER } from "./modules/gatk/haplotype_caller/haplotype_caller.nf" addParams(params.GATK_HAPLOTYPE_CALLER)
include { GATK_MARK_DUPLICATES } from "./modules/gatk/mark_duplicates/mark_duplicates.nf" addParams(params.GATK_MARK_DUPLICATES)
include { GATK_MERGE_BAM_ALIGNMENT } from "./modules/gatk/merge_bam_alignment/merge_bam_alignment.nf" addParams(params.GATK_MERGE_BAM_ALIGNMENT)
include { GATK_MERGE_VCFS } from "./modules/gatk/merge_vcfs/merge_vcfs.nf" addParams(params.GATK_MERGE_VCFS)
include { GATK_SORT_AND_FIX_TAGS } from "./modules/gatk/sort_and_fix_tags/sort_and_fix_tags.nf" addParams(params.GATK_SORT_AND_FIX_TAGS)
include { PICARD_SAM_TO_FASTQ_BWA_MEM } from "./modules/picard/sam_to_fastq_bwa_mem/sam_to_fastq_bwa_mem.nf" addParams(params.PICARD_SAM_TO_FASTQ_BWA_MEM)


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


// Prepare channel inputs for Base Recalibration

subgrouping_ch = channel.fromPath(sequence_grouping)
        .splitText()
        .map { line -> [line.tokenize(':')[0].trim(), line.trim()] }

// Prepare channel inputs for  applying BQSR

recal_scatter_counter = 0
subgrouping_unmapped_ch = channel.fromPath(sequence_grouping_unmapped)
        .splitText()
        .map { line ->
            recal_scatter_counter = recal_scatter_counter + 1
            [recal_scatter_counter, line.tokenize(':')[0].trim(), line.trim()]
        }


// Prepare channel for scattered calling intervals: input to haplotypecaller

calling_scatter_counter = 0

calling_intervals = channel.fromPath(scattered_calling_interval)
        .splitText()
        .map { line ->
            calling_scatter_counter = calling_scatter_counter + 1
            [calling_scatter_counter, line.split("/")[5], line.trim()]
        }


//================================================================================
// Define sub-workflows
//================================================================================

workflow PREPROCESSING_MAPPING {
    take:
    unmapped_bams

    main:
    PICARD_SAM_TO_FASTQ_BWA_MEM(
            unmapped_bams,
            ref_alt,
            ref_amb,
            ref_ann,
            ref_bwt,
            ref_pac,
            ref_sa,
            ref_dict,
            ref_fasta,
            ref_fasta_fai
    )

    bwa_version = BWA_GET_BWA_VERSION()

    GATK_MERGE_BAM_ALIGNMENT(
            PICARD_SAM_TO_FASTQ_BWA_MEM.out[0],
            PICARD_SAM_TO_FASTQ_BWA_MEM.out[1],
            PICARD_SAM_TO_FASTQ_BWA_MEM.out[2],
            ref_dict,
            ref_fasta,
            ref_fasta_fai,
            bwa_version
    )


    GATK_MARK_DUPLICATES(
            GATK_MERGE_BAM_ALIGNMENT.out[0],
            GATK_MERGE_BAM_ALIGNMENT.out[1]
    )

    GATK_SORT_AND_FIX_TAGS(
            GATK_MARK_DUPLICATES.out[0],
            GATK_MARK_DUPLICATES.out[1],
            ref_dict,
            ref_fasta,
            ref_fasta_fai
    )

    emit:
    GATK_SORT_AND_FIX_TAGS.out
}


//------------------

workflow QUALITY_RECALIBRATION {
    take:
    sorted_bam_and_subgroup

    main:
    GATK_BASE_RECALIBRATOR(
            sorted_bam_and_subgroup.combine(subgrouping_ch),
            ref_dict,
            ref_fasta,
            ref_fasta_fai,
            dbSNP_vcf,
            dbSNP_vcf_index,
            known_indels_mills,
            known_indels_mills_index,
            known_indels_dbSNP,
            known_indels_dbSNP_index
    )

    GATK_GATHER_BQSR_REPORTS(
            GATK_BASE_RECALIBRATOR.out.groupTuple()
    )


    GATK_APPLY_BQSR(
            sorted_bam_and_subgroup.join(GATK_GATHER_BQSR_REPORTS.out).combine(subgrouping_unmapped_ch),
            ref_dict,
            ref_fasta,
            ref_fasta_fai
    )

    GATK_GATHER_BAM_FILES(
            GATK_APPLY_BQSR.out.groupTuple()
    )

    emit:
    GATK_GATHER_BAM_FILES.out

}


//------------------

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

