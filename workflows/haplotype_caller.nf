nextflow.preview.dsl = 2

include { SAMTOOLS_CRAM_TO_BAM } from "../modules/samtools/cram_to_bam/cram_to_bam"
include { GATK_HAPLOTYPE_CALLER } from "../modules/gatk/haplotype_caller/haplotype_caller"
include { GATK_MERGE_GVCFS } from "../modules/gatk/merge_gvcfs/merge_gvcfs"



workflow {
    take:

    main:
    SAMTOOLS_CRAM_TO_BAM
    GATK_HAPLOTYPE_CALLER
    GATK_MERGE_GVCFS

    emit:

}
