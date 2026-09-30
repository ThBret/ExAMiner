process PREPARE_ALIGNMENT_INPUT {

    tag "prepare_alignment_input"
    label 'local_small'

    input:
	path alignment_source_dir

    output:
	path "prepared/*.fasta", emit: alignment_input

    script:

    def source = params.results_from_third_party ?: 'stitcher'

    if (!(source in ['stitcher', 'atram', 'patchwork'])) {
        error """
Invalid value for --results-from-third-party: '${source}'

Valid values are:
stitcher
atram
patchwork

If omitted, STITCHER results are used.
"""
    }

    """
    mkdir -p prepared

    echo "=================================================="
    echo "PREPARING ALIGNMENT INPUT"
    echo "Source: ${source}"
    echo "=================================================="

    ############################################################
    # STITCHER
    #
    # 5.3-stitcher.output/
    #   BAIT/
    #       stitcher_out/
    #           BAIT.fasta
    #
    # Each FASTA contains sample-labelled sequences.
    ############################################################

    if [ "${source}" = "stitcher" ]; then

        echo "Using STITCHER results"

        find -L ${alignment_source_dir} \\
            -type f \\
            -path "*/stitcher_out/*.fasta" |
        while read -r F; do

            BAIT=\$(basename "\$F" .fasta)
            OUT="prepared/\${BAIT}.fasta"

            echo "STITCHER: bait=\${BAIT}"
            cat "\$F" >> "\$OUT"

        done

    ############################################################
    # aTRAM
    ############################################################

    elif [ "${source}" = "atram" ]; then

        echo "Using aTRAM results"

        find -L ${alignment_source_dir} \\
            -type f \\
            -name "*.stitched_exons.fasta" |
        while read -r F; do

            BASENAME=\$(basename "\$F")
            BAIT=\$(basename "\$F" .stitched_exons.fasta)
            BAIT="\${BAIT#atram_stitcher_}"
            BAIT="\${BAIT#*.}"

            if [ -z "\$BAIT" ]; then
                echo "WARNING: Could not determine bait from \$F" >&2
                continue
            fi

            OUT="prepared/\${BAIT}.fasta"

            echo "aTRAM: bait=\${BAIT}"
            cat "\$F" >> "\$OUT"

        done

    ############################################################
    # PATCHWORK
    ############################################################

    elif [ "${source}" = "patchwork" ]; then

        echo "Using Patchwork results"

        find -L ${alignment_source_dir} \\
            -type f \\
            -path "*/nucleotide_query_sequences/*.fas" |
        while read -r F; do

            SAMPLE=\$(basename "\$(dirname "\$(dirname "\$F")")")
            BAIT=\$(basename "\$F" .fas)
            OUT="prepared/\${BAIT}.fasta"

            echo "Patchwork: sample=\${SAMPLE} bait=\${BAIT}"

            # Replace Patchwork read IDs with the sample ID
            awk -v sample="\$SAMPLE" '
                /^>/ {
                    print ">" sample
                    next
                }
                {
                    print
                }
            ' "\$F" >> "\$OUT"

        done

    fi

    ############################################################
    # VALIDATE OUTPUT
    ############################################################

    echo "=================================================="
    echo "PREPARED ALIGNMENT INPUT"
    echo "=================================================="

    FOUND=0

    for F in prepared/*.fasta; do

        [ -e "\$F" ] || continue

        BAIT=\$(basename "\$F" .fasta)
        N=\$(grep -c '^>' "\$F" || true)

        if [ "\$N" -eq 0 ]; then
            rm -f "\$F"
            continue
        fi

        echo "\${BAIT}: \${N} sequences"

        FOUND=1

    done

    if [ "\$FOUND" -eq 0 ]; then
        echo "ERROR: No valid FASTA files were produced." >&2
        exit 1
    fi

    echo "=================================================="
    echo "PREPARATION COMPLETE"
    echo "=================================================="
    """
}
