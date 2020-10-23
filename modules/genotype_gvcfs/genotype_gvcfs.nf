nextflow.enable.dsl = 2

params.output_vcf_filename
params.allow_old_rms_mapping_quality_annotation_data
params.workspace_tar

process GATK_GENOTYPE_GVCFS {
    container "us.gcr.io/broad-gatk/gatk:4.1.1.0"

    input:
    tuple path(ref_fasta), path(ref_fasta_index)
    path(ref_dict)
    path(interval)

    output:

    shell:
    '''
    set -euo pipefail

    tar -xf !{params.workspace_tar}
    WORKSPACE=$(basename !{params.workspace_tar} .tar)

    gatk --java-options -Xms8g  GenotypeGVCFs \
                                -R !{ref_fasta} \
                                -O !{output_vcf_filename} \
                                -D !{dbsnp_vcf} \
                                --only-output-calls-starting-in-intervals \
                                --use-new-qual-calculator \
                                -V gendb://$WORKSPACE \
                                -L !{interval} \
                                !{params.allow_old_rms_mapping_quality_annotation_data ? '--allow-old-rms-mapping-quality-annotation-data': '' } \
                                --merge-input-intervals
    
    
    '''
}


