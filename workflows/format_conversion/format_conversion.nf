nextflow.enable.dsl = 2

//================================================================================
// Include modules and (soft) override module-level parameters
//================================================================================

include { GATK_PAIRED_FASTQ_TO_UNMAPPED_BAM } from "../../modules/gatk/paired_fastq_to_unmapped_bam/paired_fastq_to_unmapped_bam.nf" addParams(params.GATK_PAIRED_FASTQ_TO_UNMAPPED_BAM)


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


//================================================================================
// Test workflow
//================================================================================

workflow test {
    fastq_files_list = file(params.input_fofn)

// Convert input manifest to a channel.
    fastq_params_ch = channel.fromPath(fastq_files_list)
            .splitText(keepHeader: false)
            .map { line ->
                cols = line.tokenize('\t')
                [
                        file(cols[2]), // fastq_1
                        file(cols[3]), // fastq_2
                        cols[6], // run_date
                        cols[1], // sample_name
                        cols[4], // library_name
                        cols[7], // platform_name
                        cols[5], // platform_name
                        cols[0], // readgroup_name
                        cols[8] // sequencing_center

                ]
            }


    FORMAT_CONVERSION({ fastq_params_ch })

}
