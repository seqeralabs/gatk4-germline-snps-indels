nextflow.enable.dsl = 2


pararms.recalibration_filename
params.tranches_filename
params.use_allele_specific_annotations


process GATK_INDELS_VARIANT_RECALIBRATOR {

    input:
    path(sites_only_variant_filtered_vcf)
    path(sites_only_variant_filtered_vcf_index)
    path(mills_resource_vcf)
    path(axiomPoly_resource_vcf)
    path(dbsnp_resource_vcf)
    path(mills_resource_vcf_index)
    path(axiomPoly_resource_vcf_index)
    path(dbsnp_resource_vcf_index)


    output:
    path("${params.recalibration_filename}")
    path("${params.recalibration_filename}.idx")
    path("${params.tranches_filename}")


    script:
    tranche_str = recalibration_tranche_values.collect().join(" -tranche ")
    an_str = recalibration_annotation_values.collect().join(" -an ")

    """
    set -euo pipefail

    gatk --java-options -Xms24g VariantRecalibrator \
                                -V ${sites_only_variant_filtered_vcf} \
                                -O ${params.recalibration_filename} \
                                --tranches-file ${params.tranches_filename} \
                                --trust-all-polymorphic \
                                -tranche ${tranche_str} \
                                -an ${an_str} \
                                ${params.use_allele_specific_annotations ? '--use-allele-specific-annotations' : ''} \
                                -mode INDEL \
                                --max-gaussians ${max_gaussians} \
                                -resource:mills,known=false,training=true,truth=true,prior=12 ${mills_resource_vcf} \
                                -resource:axiomPoly,known=false,training=true,truth=false,prior=10 ${axiomPoly_resource_vcf} \
                                -resource:dbsnp,known=true,training=false,truth=false,prior=2 ${dbsnp_resource_vcf}
    
    """
}
