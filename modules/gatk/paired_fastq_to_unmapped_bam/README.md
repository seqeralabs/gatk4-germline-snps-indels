# Paired fastq to unmapped bam

The pupose of module is straight-forward, it takes the `fastq` files as input and produces the `unmapped bam`. 

This particular module is imported in the `/workflow/format_conversion/format_conversion.nf` workflow and unless the parameters are overridden, the default values are used for 
`container`, `memory`, `cpus`, `java_opts` and `gatk_path`.

 
## Testing

Since DSL2 allows a module to have a worklow as well, we can test this module and optimize directives using the `test` workflow. 

This process could be tested using the `test` workflow which relies on 

- the local files in `test_data` 
- params in  `test_params.yaml`
- `test` profile in `nextflow.config`


```
nextflow run paired_fastq_to_unmapped_bam.nf -entry test -params-file test_params.yaml -profile test
```

We can use the optimum configuration without making any change to the module file, for the production workloads, by simply copy and paste the desired configs to the workflow file or workflow params.


- Workflow file e.g. `nextflow.config`

```nextflow
params.GATK_PAIRED_FASTQ_TO_UNMAPPED_BAM [
// add the optimized params
java_opts : ""
]
include { GATK_PAIRED_FASTQ_TO_UNMAPPED_BAM } from "relative_path_to/modules/gatk/paired_fastq_to_unmapped_bam/paired_fastq_to_unmapped_bam.nf" addParams (params.GATK_PAIRED_FASTQ_TO_UNMAPPED_BAM)
```

- Workflow params file e.g. `params.yaml`

```yaml

GATK_PAIRED_FASTQ_TO_UNMAPPED_BAM:
  java_opts : ""

```
