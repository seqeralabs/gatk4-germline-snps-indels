nextflow.enable.dsl = 2


//================================================================================
// Read and derive file names and location from the params
//================================================================================

ref_fasta = file(params.fasta)
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

//================================================================================
// Include modules and (soft) override module-level parameters
//================================================================================


include { GATK_APPLY_BQSR } from "../../modules/gatk/apply_bqsr/apply_bqsr.nf" addParams(params.GATK_APPLY_BQSR)
include { GATK_BASE_RECALIBRATOR } from "../../modules/gatk/base_recalibrator/base_recalibrator.nf" addParams(params.GATK_BASE_RECALIBRATOR)
include { GATK_GATHER_BAM_FILES } from "../../modules/gatk/gather_bam_files/gather_bam_files.nf" addParams(params.GATK_GATHER_BAM_FILES)
include { GATK_GATHER_BQSR_REPORTS } from "../../modules/gatk/gather_bqsr_reports/gather_bqsr_reports.nf" addParams(params.GATK_GATHER_BQSR_REPORTS)


//================================================================================
// Prepare channels
//================================================================================

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


//================================================================================
// Main workflow
//================================================================================

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


