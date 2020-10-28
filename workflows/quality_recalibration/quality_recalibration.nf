workflow QUALITY_RECALIBRATION {
    take:
    data

    main:
    baseRecalibrator(
            data.combine(subgrouping),  \
                 ref_dict,    \
                 ref_fasta,   \
                 ref_fasta_fai,   \
                 dbSNP_vcf,   \
                 dbSNP_vcf_index,     \
                known_indels_mills, \
                known_indels_mills_index,   \
                known_indels_dbSNP,     \
                known_indels_dbSNP_index   \
        )

    gatherBqsrReports(
            baseRecalibrator.out.groupTuple()
    )


    applyBQSR(
            data.join(gatherBqsrReports.out).combine(subgrouping_unmapped), \
                ref_dict,   \
                ref_fasta,  \
                ref_fasta_fai  \
        )

    gatherBamFiles(
            applyBQSR.out.groupTuple()
    )

    emit:
    gatherBamFiles.out

}


