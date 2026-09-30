process ALIGNER {
    tag "$bait"
    label 'aligning'
    publishDir "${params.outdir}/7-aligning/${params.ucode}_N${params.n_percent}S${params.s_percent}/3-alignments", mode: 'copy'

    input:
        tuple val(bait), path(exon_fastas, stageAs: 'exon_*.fasta')
        path(samples_list)

    output:
        tuple val(bait), path("${bait}_NT.fasta"), path("${bait}_AA.fasta"), optional: true, emit: alignment

    script:
    """
    mkdir -p passed

    echo "bait = ${bait}"

    echo "exon_fastas:"
    for f in ${exon_fastas}; do
        echo "  \$f"
        BLENGTH=\$(grep -v '^>' "\$f" | tr -d '\n' | wc -c)
        NCOUNT=\$(grep -v '^>' "\$f" | tr -d '\n' | tr -cd 'N' | wc -c)
        NTHRESH=\$(awk -v len=\${BLENGTH} -v pct=${params.n_percent} 'BEGIN {printf "%.0f", len*pct/100}')

        if (( NCOUNT < NTHRESH )); then
            cp \$f passed/
        fi
    done

    if [ -z "\$(ls -A passed)" ]; then
        echo "--------------------------------------------------"
        echo "LOCUS ${bait} DID NOT PASS THE MISSING DATA FILTER"
        echo "--------------------------------------------------"
        exit 0
    fi

    echo "--------------------------------------------------"
    echo "CHECKING ALIGNMENT REQUIREMENTS"
    echo "--------------------------------------------------"
    echo "passed contains:"
    ls -l passed

    cat passed/*.fasta > "${bait}.alignment_input.fasta"
    SCOUNT=\$(grep -c ">" "${bait}.alignment_input.fasta")
    TOTAL_N_SAMP=\$(grep -c '[^[:space:]]' "${samples_list}")
    STHRESH=\$(awk -v n=\${TOTAL_N_SAMP} -v pct=${params.s_percent} 'BEGIN {printf "%.0f", n*pct/100}')

    echo "--------------------------------------------------"
    echo "SAMPLE COUNT CHECK"
    echo "--------------------------------------------------"
    echo "SCOUNT (sequences in alignment input) = \${SCOUNT}"
    echo "TOTAL_N_SAMP (from samples list)      = \${TOTAL_N_SAMP}"
    echo "STHRESH (min samples required)        = \${STHRESH}"

    if (( SCOUNT < STHRESH )); then
        echo "${bait} has fewer than \${STHRESH} samples, skipping alignment" >&2
        exit 0
    fi

    mv "${bait}.alignment_input.fasta" "${bait}.fasta"
    macse -prog alignSequences -seq ${bait}.fasta
    """
}
