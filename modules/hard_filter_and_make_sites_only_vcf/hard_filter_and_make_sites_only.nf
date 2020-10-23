nextflow.enable.dsl = 2

params.variant_filtered_vcf_filename
params.excess_het_threshold
params.sites_only_vcf_filename

process GATK_HARD_FILTER_AND_MAKE_SITES_ONLY_VCF {
    container "us.gcr.io/broad-gatk/gatk:4.1.1.0"

    input:
    tuple path(vcf), path(vcf_index)
    path(interval)


    output:

    script:
    """  
    set -euo pipefail

    gatk --java-options -Xms3g VariantFiltration \
                        --filter-expression "ExcessHet > ${params.excess_het_threshold}" \
                        --filter-name ExcessHet \
                        -O ${params.variant_filtered_vcf_filename} \
                        -V ${vcf}

    gatk --java-options -Xms3g MakeSitesOnlyVcf \
                        -I ${params.variant_filtered_vcf_filename} \
                        -O ${params.sites_only_vcf_filename}
  """
}



