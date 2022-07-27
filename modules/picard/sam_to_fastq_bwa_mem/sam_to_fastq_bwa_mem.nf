params.samtools_path = "samtools"
params.bwa_path = "/usr/gitc/bwa"
params.picard_path = "/usr/gitc/picard.jar"
params.java_path = "java"
params.java_opts = "-Xms3000m"
params.compression_level = 5

process PICARD_SAM_TO_FASTQ_BWA_MEM {
    tag "${sampleId}"
    label "gitc_container"

    input:
    tuple val(sampleId), path(input_unmapped_bam)

    path(ref_alt)
    path(ref_amb)
    path(ref_ann)
    path(ref_bwt)
    path(ref_pac)
    path(ref_sa)
    path(ref_dict)
    path(ref_fasta)
    path(ref_fasta_fai)

    output:
    val(sampleId)
    path("${sampleId}.mapped.bam")
    path(input_unmapped_bam)

    script:

    """
    set -o pipefail
    set -e

	${params.java_path} -Dsamjdk.compression_level=${params.compression_level} ${params.java_opts} \
	    -jar ${params.picard_path} \
        SamToFastq \
        INPUT=${input_unmapped_bam} \
        FASTQ=/dev/out \
        INTERLEAVE=true \
        NON_PF=true \
    | \
		${params.bwa_path} mem \
		 -K 100000000 -p -v 3 -t 16 -Y ${ref_fasta} /dev/stdin -  2> >(tee ${sampleId}.bwa.stderr.log >&2) \
    | \
		${params.samtools_path} view -1 - > ${sampleId}.mapped.bam
    """


    stub:

    """
    touch ${sampleId}.mapped.bam
    """
}
