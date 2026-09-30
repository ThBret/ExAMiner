process EXONERATE_SRA {
    tag "$bait"
    label 'stitching'
    publishDir "${params.exonerate_sra_dir}/${bait}", mode: 'copy'

    input:
        tuple val(bait), path(contigs_dir), path(bait_nostop_fasta)

    output:
        tuple val(bait), path("exonerate_out"), emit: results

    script:
    """
    mkdir -p exonerate_out/${bait}

    for f in ${contigs_dir}/*_all_contigs.fasta; do
        [ -e "\$f" ] || continue

        # Skip FASTA if empty
        if ! grep -q '^>' "\$f"; then
            echo "WARNING: No sequences found in file \$f, skipping"
            continue
        fi

        S=\$(basename \$f _all_contigs.fasta)
        exonerate --useaatla yes --model protein2genome ${bait_nostop_fasta} \$f > exonerate_out/${bait}/\${S}_exonerate.out
    done
    """
}
