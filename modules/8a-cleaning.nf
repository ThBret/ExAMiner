process CLEANING {
    tag "$bait"
    label 'cleaning'
    publishDir "${params.cleaning_dir}/${params.ucode}_N${params.n_percent}S${params.s_percent}_fSTP_fFRS/8.5-filtering.input", mode: 'copy'

    input:
        tuple val(bait), path(nt_fasta)
        path outgroups_list

    output:
        tuple val(bait), path("${bait}.all.fasta"), path("${bait}.treefile"), path("${bait}.runtrees"), optional: true, emit: gene_tree
        tuple val(bait), path("${bait}.ingroup.fasta"), optional: true, emit: ingroup

    script:
    """
    export SEDA_JAVA_MEMORY=-Xmx10G

    echo "--------------------------------------------------"
    echo "PREPARING LOCUS ${bait}"
    echo "--------------------------------------------------"

    if grep -q '!' ${nt_fasta}; then
        sed 's/!/?/g' ${nt_fasta} > input.fasta
    else
        cp ${nt_fasta} input.fasta
    fi

    ##################################################
    ### CHECK INTERNAL STOP CODONS
    ##################################################

    echo "--------------------------------------------------"
    echo "CHECKING FOR INTERNAL STOP CODONS AND FRAMESHIFTS"
    echo "--------------------------------------------------"

    echo "Running SEDA"

    mkdir -p stopfiltered/with.frameshift stopfiltered/without.frameshift

    set +e
    seda filtering -if input.fasta -od stopfiltered --remove-with-in-frame-stop-codons --in-disk-processing
    SEDA_EXIT=\$?
    set -e

    PREFILT=\$(grep -c ">" input.fasta)

    # SEDA sometimes exits successfully but produces no output when all sequences
    # are removed by the in-frame stop codon filter.
    if [[ \$SEDA_EXIT -ne 0 ]] || [[ ! -s stopfiltered/input.fasta ]]; then
        echo "Locus ${bait}: no sequences survived stop-codon filtering" >&2
        exit 0
    fi

    POSTFILT=\$(grep -c ">" stopfiltered/input.fasta)

    DIFF=\$((PREFILT - POSTFILT))

    if (( DIFF >= ${params.stop_thresh} )); then
        echo "Locus ${bait} internal stop codon filter: FAILED" >&2

        if grep -q '?' stopfiltered/input.fasta; then
            mv stopfiltered/input.fasta stopfiltered/with.frameshift/${bait}_NT.fasta
            sed -i 's/?/!/g' stopfiltered/with.frameshift/${bait}_NT.fasta
        else
            mv stopfiltered/input.fasta stopfiltered/without.frameshift/${bait}_NT.fasta
        fi

        exit 0
    else
        echo "Locus ${bait} internal stop codon filter: PASSED"

    fi

    ##################################################
    ### FRAMESHIFT CHECK
    ##################################################

    cp stopfiltered/input.fasta frs_check.fasta

    if grep -q '?' frs_check.fasta; then
        echo "Locus ${bait} frameshift filter: FAILED"
        echo "Check for erroneous frameshift at locus start/end positions due to missing data"
        sed -i 's/?/!/g' frs_check.fasta
        exit 0
    else
        echo "Locus ${bait} frameshift filter: PASSED"
    fi

    ##################################################
    ### CLIPKIT TRIMMING
    ##################################################

    echo "--------------------------------------------------"
    echo "TRIMMING WITH CLIPKIT"
    echo "--------------------------------------------------"

    mkdir -p trimlogs

    clipkit frs_check.fasta \\
        --output ${bait}_NT.trim.fasta \\
        --mode smart-gap \\
        -s nt \\
        --codon \\
        --log

    mv frs_check.fasta ${bait}_NT.fasta
    mv ${bait}_NT.trim.fasta.log trimlogs/

    ##################################################
    ### CHECK OUTGROUPS
    ##################################################

    awk '
        BEGIN{
            while((getline < "'"${outgroups_list}"'")>0) wanted[\$1]=1
        }
        /^>/{
            id=substr(\$0,2)
            keep = wanted[id]
            next
        }
        keep && /[ACGT]/{
            print id
            delete wanted[id]
            keep=0
        }
    ' ${bait}_NT.trim.fasta > tmp.list.txt
    
    [[ -s tmp.list.txt ]] || {
        echo "Locus ${bait} - no valid outgroup with nucleotides after trimming, skipping gene tree" >&2
        exit 0
    }

    OUTGROUP_CSV=\$(paste -sd, tmp.list.txt)
    OUTGROUP_SPACE=\$(paste -sd' ' tmp.list.txt)

    ##################################################
    ### IQTREE GENE TREE
    ##################################################

    echo "--------------------------------------------------"
    echo "RUNNING IQTREE FOR LOCUS ${bait}"
    echo "--------------------------------------------------"

    # Create partition file
    SEQLEN=\$(seqkit stats -T ${bait}_NT.trim.fasta | tail -n 1 | awk '{print \$8}')

    printf '%s\n' \
        "DNA, Subset1 = 1-\${SEQLEN}\\3" \
        "DNA, Subset2 = 2-\${SEQLEN}\\3" \
        "DNA, Subset3 = 3-\${SEQLEN}\\3" \
        > ${bait}_NT.trim.prt

    # Run iqtree
    iqtree \\
        -s ${bait}_NT.trim.fasta \\
        -p ${bait}_NT.trim.prt \\
        -m MFP \\
        --merge \\
        --seqtype DNA \\
        -o \${OUTGROUP_CSV} \\
        --runs ${params.iqtree_runs} \\
        -B ${params.iqtree_boot} \\
        --nstop 500 \\
        --bnni \\
        -T ${task.cpus}

    ##################################################
    ### SAVE IQTREE OUTPUTS
    ##################################################

    cp ${bait}_NT.trim.fasta ${bait}.all.fasta
    cp ${bait}_NT.trim.prt.treefile ${bait}.treefile
    cp ${bait}_NT.trim.prt.runtrees ${bait}.runtrees

    ##################################################
    ### REMOVE OUTGROUPS
    ##################################################

    AMAS.py remove \\
        -i ${bait}.all.fasta \\
        -f fasta \\
        -d dna \\
        -x \${OUTGROUP_SPACE}

    mv reduced_${bait}.all.fasta-out.fas ${bait}.ingroup.fasta
    """
}
