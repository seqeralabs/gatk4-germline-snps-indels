nextflow.enable.dsl = 2


params.scattered
params.output_base_name

process PICARD_CROSS_CHECK_FINGERPRINT {
    container "us.gcr.io/broad-gotc-prod/gatk4-joint-genotyping:yf_fire_crosscheck_picard_with_nio_fast_fail_fast_sample_map"

    input:
    path(gvcf_paths)
    path(vcf_paths)
    path(sample_name_map)
    path(haplotype_database)
    path(expected_inconclusive_samples)


    output:
    path("${output_name}")


    shell:
    output_name = output_base_name + ".fingerprintcheck"
    expected_inconclusive_samples = {sep='" "' expected_inconclusive_samples}

    '''
    set -eu

    gvcfInputsList=!{write_lines(gvcf_paths)}
    vcfInputsList=!{write_lines(vcf_paths)}

    cp $gvcfInputsList gvcf_inputs.list
    cp $vcfInputsList vcf_inputs.list

    java -Dpicard.useLegacyParser=false -Xms!{memMb - 512}m \
      -jar /usr/gitc/PicardPublicWithCrosscheckNIOandSampleMapping.jar \
      CrosscheckFingerprints \
      --INPUT gvcf_inputs.list \
      --SECOND_INPUT vcf_inputs.list \
      --HAPLOTYPE_MAP !{haplotype_database} \
      --INPUT_SAMPLE_FILE_MAP !{sample_name_map} \
      --CROSSCHECK_BY SAMPLE \
      --CROSSCHECK_MODE CHECK_SAME_SAMPLE \
      --NUM_THREADS !{cpu} \
      --SKIP_INPUT_READABLITY_TEST \
      !{scattered ? '--EXIT_CODE_WHEN_MISMATCH 0': '' } \
      --OUTPUT !{output_name}

    if !{scattered}; then
        # UNEXPECTED_MATCH is not possible with CHECK_SAME_SAMPLE
        matches=$(grep "EXPECTED_MATCH" !{output_name} | wc -l)

        # check inconclusive samples
        expectedInconclusiveSamples=("!{expected_inconclusive_samples}")
        inconclusiveSamplesCount=0
        inconclusiveSamples=($(grep 'INCONCLUSIVE' !{output_name} | cut -f 1))
        
        for sample in ${inconclusiveSamples[@]}; do
            if printf '%s\n' ${expectedInconclusiveSamples[@]} | grep -P '^'${sample}'$'; then
            inconclusiveSamplesCount=$((inconclusiveSamplesCount+1))
            fi
        done

        total_matches=$((inconclusiveSamplesCount + matches))
        
        if [[ ${total_matches} -eq !{num_gvcfs} ]]; then
            >&2 echo "Found the correct number of matches (!{num_gvcfs}) for this shard"
        else
            >&2 echo "ERROR: Found $total_matches 'EXPECTED_MATCH' records, but expected !{num_gvcfs}"
            exit 1
        fi
    fi
    '''
}
