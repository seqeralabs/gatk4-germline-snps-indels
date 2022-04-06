nextflow.enable.dsl = 2


//================================================================================
// Read and derive file names and location from the params
//================================================================================

ref_fasta = file(params.fasta)
ref_alt = file("${params.fasta}.64.alt")
ref_amb = file("${params.fasta}.64.amb")
ref_ann = file("${params.fasta}.64.ann")
ref_bwt = file("${params.fasta}.64.bwt")
ref_pac = file("${params.fasta}.64.pac")
ref_sa = file("${params.fasta}.64.sa")
ref_fasta_fai = file("${params.fasta}.fai")
ref_dict = file(params.fasta.replace(".fasta", ".dict"))

//================================================================================
// Include modules and (soft) override module-level parameters
//================================================================================

include { PICARD_SAM_TO_FASTQ_BWA_MEM } from "../../modules/picard/sam_to_fastq_bwa_mem/sam_to_fastq_bwa_mem.nf" addParams(params.PICARD_SAM_TO_FASTQ_BWA_MEM)
include { BWA_GET_BWA_VERSION } from "../../modules/bwa/get_bwa_version/get_bwa_version.nf"
include { GATK_MERGE_BAM_ALIGNMENT } from "../../modules/gatk/merge_bam_alignment/merge_bam_alignment.nf" addParams(params.GATK_MERGE_BAM_ALIGNMENT)
include { GATK_MARK_DUPLICATES } from "../../modules/gatk/mark_duplicates/mark_duplicates.nf" addParams(params.GATK_MARK_DUPLICATES)
include { GATK_SORT_AND_FIX_TAGS } from "../../modules/gatk/sort_and_fix_tags/sort_and_fix_tags.nf" addParams(params.GATK_SORT_AND_FIX_TAGS)


//================================================================================
// Main workflow
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

