workflow preprocessing_mapping {
    take: data

    main:
    samToFastqBwaMem(data,ref_alt,ref_amb,ref_ann,ref_bwt,ref_pac,ref_sa,ref_dict,ref_fasta,ref_fasta_fai)

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

    markDuplicates(mergeBamAlignment.out[0],mergeBamAlignment.out[1])

    sortAndFixTags(markDuplicates.out[0],markDuplicates.out[1],ref_dict,ref_fasta,ref_fasta_fai)

    emit:
    sortAndFixTags.out
}


