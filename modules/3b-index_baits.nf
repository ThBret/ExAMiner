process INDEX_BAITS_STOP {
    label 'slurm_small'
    publishDir "${params.infodir}", mode: 'copy'

    input:
        path baits_stop_dir

    output:
        tuple path("baits_stop_combined.fasta"),
              path("baits_stop_combined.fasta.fai"),
              path("baits_stop_id_map.tsv"), emit: baits_stop_indexed

    script:
    """
    : > baits_stop_id_map.tsv
    : > baits_stop_combined.fasta

    for f in ${baits_stop_dir}/*.fasta; do
        BAIT=\$(basename "\$f" .fasta | sed -E 's/^.*_//')
        HEADER=\$(grep -m1 '^>' "\$f" | sed 's/^>//' | awk '{print \$1}')
        echo -e "\${BAIT}\\t\${HEADER}" >> baits_stop_id_map.tsv
        cat "\$f" >> baits_stop_combined.fasta
    done

    samtools faidx baits_stop_combined.fasta
    """
}
