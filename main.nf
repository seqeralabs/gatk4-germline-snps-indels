nextflow.enable.dsl = 2

//--------------


java_opt_samtofastq = "-Xms3000m"
java_opt_mergebams = "-Xms6000m"
java_opt_markdups = "-Xms4000m"
java_opt_sort = "-Xms4000m"
java_opt_fix = "-Xms500m"
java_opt_baserecal = "-Xms4000m"
java_opt_bqsrreport = "-Xms3000m"
java_opt_applybqsr = "-Xms3000m"
java_opt_gatherbams = "-Xms2000m"
java_opt_haplotype = "-Xmx4G"
java_opt_mergevcfs = "-Xmx4G"
java_opt_genomicsDBimport = "-Xmx4G"


compression_level = 5


unmapped_bams_list = "s3://nf-work-bucket/unmapped_bams_3.tsv"
outdir = "s3://nf-work-bucket/results/M5-6-samples"

fasta = "s3://broad-references/hg38/v0/Homo_sapiens_assembly38.fasta"
dbSNP_vcf = "s3://broad-references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf"
known_indels_mills = "s3://broad-references/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
known_indels_dbSNP = "s3://broad-references/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz"
sequence_grouping = "s3://da-gatk-nf-dev/sequence_grouping.txt"
sequence_grouping_unmapped = "s3://da-gatk-nf-dev/sequence_grouping_with_unmapped.txt"
scattered_calling_interval = "s3://da-gatk-nf-dev/hg38_wgs_scattered_calling_intervals.txt"

//--------------

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

include { GATK_APPLY_BQSR } from "../modules/gatk/apply_bqsr/apply_bqsr.nf"
include { GATK_BASE_RECALIBRATOR } from "../modules/gatk/base_recalibrator/base_recalibrator.nf"
include { GATK_GATHER_BAM_FILES } from "../modules/gatk/gather_bam_files/gather_bam_files.nf"
include { GATK_GATHER_BQSR_REPORTS } from "../modules/gatk/gather_bqsr_reports/gather_bqsr_reports.nf"
include { GATK_HAPLOTYPE_CALLER } from "../modules/gatk/haplotype_caller/haplotype_caller.nf"
include { GATK_MARK_DUPLICATES } from "../modules/gatk/mark_duplicates/mark_duplicates.nf"
include { GATK_MERGE_BAM_ALIGNMENT } from "../modules/gatk/merge_bam_alignment/merge_bam_alignment.nf"
include { GATK_MERGE_VCFS } from "../modules/gatk/merge_vcfs/merge_vcfs.nf"
include { GATK_SORT_AND_FIX_TAGS } from "../modules/gatk/sort_and_fix_tags/sort_and_fix_tags.nf"
include { PICARD_SAM_TO_FASTQ_BWA_MEM } from "../modules/picard/sam_to_fastq_bwa_mem/sam_to_fastq_bwa_mem.nf"
include { UTILS_CREATE_SEQUENCE_GROUPING } from "../modules/utils/create_sequence_grouping/create_sequence_grouping.nf"
include { UTILS_GET_BWA_VERSION } from "../modules/utils/get_bwa_version/get_bwa_version.nf"


/*

If we use a tab-delimited files with sampleId in one column and path to unaligned bam in the second,
it will eliminate the need make assumptions about base filename structure

*/
unmapped_bams_channel = channel.fromPath(unmapped_bams)
        .splitText()
        .map { line -> [line.tokenize("\t")[0], file(line.tokenize("\t")[1].trim())] }


// Prepare channel inputs for Base Recalibration

subgrouping = channel.fromPath(sequence_grouping)
        .splitText()
        .map { line -> [line.tokenize(':')[0].trim(), line.trim()] }

// Prepare channel inputs for  applying BQSR

recal_scatter_counter = 0
subgrouping_unmapped = channel.fromPath(sequence_grouping_unmapped)
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
            [calling_scatter_counter, line.split("/")[6], line.trim()]
        }


/*
*
################################################ EXECUTION #############################################################
*
*/

workflow preprocessing_mapping {
    take:
    data

    main:
    samToFastqBwaMem(data, ref_alt, ref_amb, ref_ann, ref_bwt, ref_pac, ref_sa, ref_dict, ref_fasta, ref_fasta_fai)

    bwa_version = getBwaVersion()

    mergeBamAlignment(
            samToFastqBwaMem.out[0],
            samToFastqBwaMem.out[1],
            samToFastqBwaMem.out[2],
            ref_dict,
            ref_fasta,
            ref_fasta_fai,
            bwa_version
    )

    markDuplicates(mergeBamAlignment.out[0], mergeBamAlignment.out[1])

    sortAndFixTags(markDuplicates.out[0], markDuplicates.out[1], ref_dict, ref_fasta, ref_fasta_fai)

    emit:
    sortAndFixTags.out
}


workflow quality_recalibraiton {
    take:
    data

    main:
    baseRecalibrator(
            data.combine(subgrouping),         \
                        ref_dict,           \
                        ref_fasta,          \
                        ref_fasta_fai,          \
                        dbSNP_vcf,          \
                        dbSNP_vcf_index,            \
                        known_indels_mills,         \
                        known_indels_mills_index,           \
                        known_indels_dbSNP,             \
                        known_indels_dbSNP_index           \
        )

    gatherBqsrReports(
            baseRecalibrator.out.groupTuple()
    )


    applyBQSR(
            data.join(gatherBqsrReports.out).combine(subgrouping_unmapped),         \
                        ref_dict,           \
                        ref_fasta,          \
                        ref_fasta_fai          \
        )

    gatherBamFiles(
            applyBQSR.out.groupTuple()
    )

    emit:
    gatherBamFiles.out

}


workflow variant_discovery {
    take:
    data

    main:

    haplotypeCaller(
            data.combine(calling_intervals),         \
                        ref_dict,           \
                        ref_fasta,          \
                        ref_fasta_fai         \
        )

    mergeVCFs(
            haplotypeCaller.out.groupTuple()
    )
    emit:
    mergeVCFs.out


}


workflow {

    preprocessing_mapping(unmapped_bams_channel)

    quality_recalibraiton(preprocessing_mapping.out)

    variant_discovery(quality_recalibraiton.out)

    variant_discovery.out.map {
        sampleId, vcfFile ->
            "${sampleId}\ts3://${vcfFile}"
    }
            .collectFile(
                    name: 'merged_vcfs.tsv', newLine: true, storeDir: "${params.outdir}"
            )
}

