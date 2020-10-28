workflow variant_discovery {
    take: data

    main:

    haplotypeCaller(
            data.combine(calling_intervals), \
                ref_dict,   \
                ref_fasta,  \
                ref_fasta_fai \
        )

    mergeVCFs(
            haplotypeCaller.out.groupTuple()
    )
    emit:
    mergeVCFs.out


}

