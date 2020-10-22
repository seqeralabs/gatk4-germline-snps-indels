import java.nio.file.Paths

nextflow.enable.dsl = 2

include { SAMTOOLS_CRAM_TO_BAM } from "../modules/samtools/cram_to_bam/cram_to_bam"
include { GATK_HAPLOTYPE_CALLER } from "../modules/gatk/haplotype_caller/haplotype_caller"
include { GATK_MERGE_GVCFS } from "../modules/gatk/merge_gvcfs/merge_gvcfs"


make_gvcf = true
make_bamout = true

def is_cram(input_bam_path) {
    return file(input_bam_path).getExtension() == "cram"
}

def get_sample_basename(input_bam_path) {
    if (is_cram(input_bam_path)) {
        sample_basename = file(input_bam_path).getBaseName() + ".cram"
    } else {
        sample_basename = file(input_bam_path).getBaseName() + ".bam"
    }

    return sample_basename
}

def vcf_basename = get_sample_basename

def output_suffix = make_gvcf ? ".g.vcf.gz" : ".vcf.gz"

def output_filename = vcf_basename + output_suffix

workflow {
    take:

    ref_fasta_ch = Channel.value([Paths.get("${baseDir}/../test_data/Homo_sapiens_assembly38.fasta"),
                                  Paths.get("${baseDir}/../test_data/Homo_sapiens_assembly38.fasta.fai")])

    ref_dict_ch = Channel.value(Paths.get("${baseDir}/../test_data/Homo_sapiens_assembly38.dict"))


    input_bam_ch = Channel.fromPath(["${baseDir}/../test_data/*.bam",
                                     "${baseDir}/../test_data/*.bai"]).buffer(size: 2)

    main:
    SAMTOOLS_CRAM_TO_BAM(
            ref_fasta_ch,
            ref_dict_ch,
            input_cram_ch
    )

//    GATK_HAPLOTYPE_CALLER
//    GATK_MERGE_GVCFS

}
