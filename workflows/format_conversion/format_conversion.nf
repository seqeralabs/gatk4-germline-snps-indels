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

