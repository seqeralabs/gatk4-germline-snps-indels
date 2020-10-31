nextflow.enable.dsl = 2

params.gitc_path = "/usr/gitc"
params.java_opts = "-Xms3000m"
params.compression_level = 5

process PICARD_SAM_TO_FASTQ_BWA_MEM {
    tag "${sampleId}"

    container "broadinstitute/genomes-in-the-cloud:2.3.1-1512499786"
    memory "16 GB"
    cpus 16

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
    path "*.mapped.bam"
    path(input_unmapped_bam)

    script:

    """
    set -o pipefail
    set -e

	java -Dsamjdk.compression_level=${params.compression_level} ${params.java_opts} \
	    -jar ${params.gitc_path}/picard.jar \
        SamToFastq \
        INPUT=${input_unmapped_bam} \
        FASTQ=/dev/stdout \
        INTERLEAVE=true \
        NON_PF=true \
    | \
		${params.gitc_path}/bwa mem \
		 -K 100000000 -p -v 3 -t 16 -Y ${ref_fasta} /dev/stdin -  2> >(tee ${sampleId}.bwa.stderr.log >&2) \
    | \
		samtools view -1 - > ${sampleId}.mapped.bam
    """
}
