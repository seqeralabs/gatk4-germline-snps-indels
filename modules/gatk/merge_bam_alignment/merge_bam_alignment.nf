process mergeBamAlignment {
    tag "${sampleId}"

    memory '16 GB'
    cpus 16

    container "broadinstitute/gatk:4.1.8.1"

    input:

    val(sampleId)
    path(input_mapped_bam)
    path(input_unmapped_bam)

    path(ref_dict)
    path(ref_fasta)
    path(ref_fasta_fai)

    val(bwa_version)

    output:
    val(sampleId)
    path "*.mapped.merged.bam"

    script:
    // DRY: Move this to config or params file
    gatk_path = "/gatk/gatk"
    bwa_commandline = "bwa mem -K 100000000 -p -v 3 -t 8 -Y ${ref_fasta}"

    """
    ${gatk_path} --java-options "-Dsamjdk.compression_level=${params.compression_level} ${params.java_opt_mergebams}" \
    MergeBamAlignment  \
    --VALIDATION_STRINGENCY SILENT \
    --EXPECTED_ORIENTATIONS FR \
    --ATTRIBUTES_TO_RETAIN X0 \
    --ALIGNED_BAM ${input_mapped_bam} \
    --UNMAPPED_BAM ${input_unmapped_bam} \
    --OUTPUT ${sampleId}.mapped.merged.bam \
    --REFERENCE_SEQUENCE ${ref_fasta} \
    --PAIRED_RUN true \
    --SORT_ORDER "unsorted" \
    --IS_BISULFITE_SEQUENCE false \
    --ALIGNED_READS_ONLY false \
    --CLIP_ADAPTERS false \
    --MAX_RECORDS_IN_RAM 2000000 \
    --ADD_MATE_CIGAR true \
    --MAX_INSERTIONS_OR_DELETIONS -1 \
    --PRIMARY_ALIGNMENT_STRATEGY MostDistant \
    --PROGRAM_RECORD_ID "bwamem" \
    --PROGRAM_GROUP_VERSION "${bwa_version}" \
    --PROGRAM_GROUP_COMMAND_LINE "${bwa_commandline}" \
    --PROGRAM_GROUP_NAME "bwamem" \
    --UNMAPPED_READ_STRATEGY COPY_TO_TAG \
    --ALIGNER_PROPER_PAIR_FLAGS true \
    --UNMAP_CONTAMINANT_READS true
    """
}