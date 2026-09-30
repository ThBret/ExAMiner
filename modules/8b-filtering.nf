process FILTERING {
    label 'filtering'
    publishDir "${params.final_out_dir}", mode: 'copy'

    input:
        path all_gene_trees   // collect()'d list of *.all.fasta / *.ingroup.fasta / *.treefile / *.runtrees from CLEANING
        path physsu_script

    output:
        path "alignments/*.all.fasta", emit: alignments_all
        path "alignments/*.ingroup.fasta", emit: alignments_ingroup
        path "trees/*.ML.treefile", emit: ml_trees
        path "stats/physsur/whitelist.txt", emit: whitelist
        path "stats/postfilter_stats_all.txt", emit: postfilter_stats
        path "stats/*", emit: stats

    script:
    """
    mkdir -p stats alignments trees

    # Get general alignment statistics - prefiltering
    AMAS.py summary -i *.all.fasta -f fasta -d dna -o stats/prefilter_stats_all.txt -c ${task.cpus}
    AMAS.py summary -i *.ingroup.fasta -f fasta -d dna -o stats/prefilter_stats_ingroup.txt -c ${task.cpus}

    # Filter using PHYSSUR
    Rscript ${physsu_script} \\
        -i . \\
        -o stats/physsur \\
        -r \\
        -b ${params.phys_boot} \\
        -p ${params.phys_prop}

    # Get general alignment statistics - postfiltering
    for b in \$(cat stats/physsur/whitelist.txt); do
        cp \${b}.all.fasta alignments/
        cp \${b}.ingroup.fasta alignments/
        cp \${b}.treefile trees/\${b}.ML.treefile
    done

    AMAS.py summary -i alignments/*.all.fasta -f fasta -d dna -o stats/postfilter_stats_all.txt -c ${task.cpus}
    AMAS.py summary -i alignments/*.ingroup.fasta -f fasta -d dna -o stats/postfilter_stats_ingroup.txt -c ${task.cpus}
    """
}
