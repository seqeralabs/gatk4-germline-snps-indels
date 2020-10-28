// From the main.nf unnamed workflow

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
