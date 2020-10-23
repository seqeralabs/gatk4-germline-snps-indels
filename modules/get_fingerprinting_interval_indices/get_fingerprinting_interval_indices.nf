nextflow.enable.dsl = 2

process GATK_GET_FINGERPRINTING_INTERVAL_INDICES {
    container "us.gcr.io/broad-gatk/gatk:4.1.1.0"

    input:
    path(unpadded_intervals)
    path(haplotype_database)

    output:
    path("indices.out")
    path("all.sorted.interval_list")
    path("all.interval_list")
    path("hdb.interval_list")

    shell:
    ''' 
    set -xeo pipefail

    function rename_intervals(){
        interval_list=$1
        name=$2

        awk 'BEGIN{FS=IFS="\t";OFS="\t";} $0~"^@"{print;next;} $0~"#CHROM"{next;} {$5="'$name'"; print}' $interval_list
    }
    
    export -f rename_intervals

    function hdb_to_interval_list(){
        input=$1

        awk 'BEGIN{IFS="\t";OFS="\t";} $0~"^@"{print;next;} $0~"#CHROM"{next;} {print $1,$2,$2,"+","interval-"NR}' $1
    }

    function rename_scatter(){
        file=$1
        number=$(echo $file | sed -E 's|([0-9]+)-scattered\.interval.*|\1|')
        rename_intervals $file $number > scattered.renamed.$number.interval_list
    }
    export -f rename_scatter

    # rename the intervals within each interval_list according to the number in the name of the list

    cp ~{sep=' ' unpadded_intervals} ./

    cat ~{write_lines(unpadded_intervals)} | xargs -n1 basename | xargs -I{} bash -c 'rename_scatter $@' _ {}

    #find the first header
    find . -name "scattered.renamed.*.interval_list" | head -n1 | xargs cat | grep '^@' > all.interval_list

    # concatenate the resulting intervals (with no header)
    find . -name "scattered.renamed.*.interval_list"  | xargs cat | grep -v '^@' >> all.interval_list

    # convert the Haplotype_database to an interval_list
    hdb_to_interval_list ~{haplotype_database} > hdb.interval_list

    # find the intervals that overlap the haplotype_database
    gatk IntervalListTools \
      -ACTION OVERLAPS \
      -O all.sorted.interval_list \
      -I all.interval_list \
      -SI hdb.interval_list

    if grep -v '^@' all.sorted.interval_list; then
      grep -v '^@' all.sorted.interval_list | awk '{FS="\t"; print $5}' | uniq > indices.out
    else
      touch indices.out
    fi
  '''

}
