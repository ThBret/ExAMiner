process PARTITIONFINDER {
    tag "$bait"
    label 'partitionfinder'
    publishDir "${params.partfi_dir}/${bait}", mode: 'copy'

    input:
        tuple val(bait), path(alignment_fasta), path(postfilter_stats), path(template_cfg)

    output:
        tuple val(bait), path("analysis/best_scheme.txt"), optional: true, emit: best_scheme
        path "${bait}.all.phy", optional: true, emit: phylip
        path "partition_finder.cfg", optional: true, emit: cfg

    script:
    """
    echo "=================================================="
    echo "RUNNING MODEL SELECTION FOR LOCUS ${bait}"
    echo "=================================================="

    ##################################################
    ### CONVERT ALIGNMENT TO RELAXED PHYLIP
    ##################################################

    goalign reformat phylip \\
        -i ${alignment_fasta} \\
        --no-block --one-line \\
        -o ${bait}.all.phy

    ##################################################
    ### LOOK UP ALIGNMENT LENGTH FROM AMAS STATS
    ##################################################

    LEN=\$(grep "${bait}" ${postfilter_stats} | awk '{print \$3}')

    if [ -z "\${LEN}" ]; then
        echo "WARNING: no alignment length found for ${bait} in ${postfilter_stats}, skipping" >&2
        exit 0
    fi

    echo "Alignment length for ${bait}: \${LEN}"

    ##################################################
    ### BUILD PARTITIONFINDER CONFIG FROM TEMPLATE
    ##################################################

    sed \\
        -e "s/REPLACEN/${bait}.all.phy/g" \\
        -e "s/REPLACEL/\${LEN}/g" \\
        ${template_cfg} > partition_finder.cfg

    ##################################################
    ### RUN PARTITIONFINDER
    ##################################################

    set +e
    partitionfinder -v --processes ${task.cpus} .
    STATUS=\$?
    set -e

    if [ \${STATUS} -ne 0 ] || [ ! -s "analysis/best_scheme.txt" ]; then
        echo "WARNING: PartitionFinder did not produce a best_scheme.txt for ${bait} (exit=\${STATUS})" >&2
        exit 0
    fi

    echo "Finished ${bait} (PartitionFinder exit=\${STATUS})"
    """
}
