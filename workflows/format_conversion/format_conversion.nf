#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

//================================================================================
// Include modules and (soft) override module-level parameters
//================================================================================

include { GATK_PAIRED_FASTQ_TO_UNMAPPED_BAM } from "../../modules/gatk/paired_fastq_to_unmapped_bam/paired_fastq_to_unmapped_bam.nf"


//================================================================================
// Prepare channels
//================================================================================

fastq_files_list = file(params.input_fofn)

// Convert input manifest to a channel.
fastq_params_ch = channel.fromPath(fastq_files_list)
        .splitText(keepHeader: false)
        .map { line ->
            cols = line.tokenize('\t')
            [
                    file(cols[2]),
                    file(cols[3]),
                    cols[6],
                    cols[1],
                    cols[4],
                    cols[7],
                    cols[5],
                    cols[0],
                    cols[8]

            ]
        }


//================================================================================
// Define sub-workflows
//================================================================================

workflow FORMAT_CONVERSION {
    take:
    fastq_files

    main:
    GATK_PAIRED_FASTQ_TO_UNMAPPED_BAM(fastq_files)
    emit:
    GATK_PAIRED_FASTQ_TO_UNMAPPED_BAM.out
}

//================================================================================
// Main workflow
//================================================================================

workflow {
    FORMAT_CONVERSION({ fastq_params_ch })
    FORMAT_CONVERSION.out
            .map { sample_name, unmapped_bam ->
                "${sample_name}\ts3://${unmapped_bam}"
            }
            .collectFile(
                    name: 'unmapped_bams.tsv', newLine: true, storeDir: "${params.outdir}"
            )
}


