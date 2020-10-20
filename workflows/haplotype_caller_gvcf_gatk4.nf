nextflow.preview.dsl = 2

//===========================
// Module imports
//===========================

include { SAMTOOLS_CRAM_TO_BAM } from "../modules/samtools/cram_to_bam/cram_to_bam"
include { GATK_HAPLOTYPE_CALLER } from "../modules/gatk/haplotype_caller/haplotype_caller"
include { GATK_MERGE_GVCFS } from "../modules/gatk/merge_gvcfs/merge_gvcfs"


//===========================
// Workflow parameters
//===========================

//===========================
// Workflow definition
//===========================

workflow joint_genotyping {
    take:

    main:

}
