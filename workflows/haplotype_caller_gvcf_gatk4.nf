nextflow.preview.dsl = 2

// Import modules
include { CRAM_TO_BAM } from "../modules/cram_to_bam/cram_to_bam"
include { HAPLOTYPE_CALLER } from "../modules/haplotype_caller/haplotype_caller"
include { MERGE_GVCFS } from "../merge_gvcfs/merge_gvcfs"


// Parameters


// Workflow
workflow joint_genotyping {
    take:

    main:

}
