nextflow.enable.dsl = 2

params.output_vcf_filename


process GATK_GNARLY_GENOTYPER {
    container "us.gcr.io/broad-gatk/gatk:4.1.1.0"

    input:
    tuple path(ref_fasta), path(ref_fasta_index)
    path(ref_dict)
    path(interval)

    output:
    tuple path("${params.output_vcf_filename}"), path("${params.output_vcf_filename}.tbi")
    path("annotationDB.vcf.gz")

    script:
    '''
    set -e

    tar -xf ~{workspace_tar}
    WORKSPACE=$( basename ~{workspace_tar} .tar)

    # use a query.json to set some params that aren't exposed -- ewwwww
    cat <<EOF > $WORKSPACE/query.json
    {
        "scan_full": true,
        "workspace": "genomicsdb",
        "array": "genomicsdb_array",
        "vid_mapping_file": "genomicsdb/vidmap.json",
        "callset_mapping_file": "genomicsdb/callset.json",
        "reference_genome": "/cromwell_root/broad-references/hg38/v0/Homo_sapiens_assembly38.fasta",
        "max_diploid_alt_alleles_that_can_be_genotyped": 6,
        "produce_GT_field": true
    }
    EOF

    gatk --java-options -Xms8g GnarlyGenotyper \
                                -R !{ref_fasta} \
                                -O !{output_vcf_filename} \
                                --output-database-name annotationDB.vcf.gz \
                                -D !{dbsnp_vcf} \
                                --only-output-calls-starting-in-intervals \
                                --use-new-qual-calculator \
                                -V gendb://$WORKSPACE \
                                -L !{interval} \
                                -stand-call-conf 10 \
                                --merge-input-intervals
    '''
}


