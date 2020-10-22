nextflow.enable.dsl = 2

include { SAMTOOLS_CRAM_TO_BAM } from "../modules/samtools/cram_to_bam/cram_to_bam"
include { GATK_HAPLOTYPE_CALLER } from "../modules/gatk/haplotype_caller/haplotype_caller"
include { GATK_MERGE_GVCFS } from "../modules/gatk/merge_gvcfs/merge_gvcfs"


params.deeptools_bamcompare_opts = ["--binsize 50"]
params.grch38_blood_vs_cfDNA = "$baseDir/test_data/features/GRCh38.blood.vs.cfDNA.v1.bed"


bam_readcount_I_features_ch = Channel.fromPath([params.grch38_blood_vs_cfDNA,
                                                params.H110TM180403])


workflow {
    take:
    bam_readcount_I_features_ch

    main:
    SAMTOOLS_CRAM_TO_BAM
    GATK_HAPLOTYPE_CALLER
    GATK_MERGE_GVCFS

}
