process MRBAYES {
    tag "$bait"
    label 'mrbayes'
    publishDir "${params.mrbayes_dir}/${bait}", mode: 'copy', pattern: "*.{nex,log}"
    publishDir "${params.bi_trees_dir}", mode: 'copy', pattern: "*.BI.treefile"

    input:
        tuple val(bait), path(alignment_fasta), path(best_scheme)
        path outgroups_list

    output:
        path "${bait}.BI.treefile", optional: true, emit: bi_tree
        path "${bait}.all.nex", optional: true, emit: nexus
        path "${bait}.log", optional: true, emit: log

    script:
    """
    echo "=================================================="
    echo "PREPARING MRBAYES INPUT FOR LOCUS ${bait}"
    echo "=================================================="

    ##################################################
    ### CONVERT ALIGNMENT TO NEXUS + APPEND MRBAYES BLOCK
    ##################################################

    goalign reformat nexus -i ${alignment_fasta} -o ${bait}.all.nex
    sed -i 's/format datatype=dna/format missing=${params.examiner_missing} gap=- datatype=dna/g' ${bait}.all.nex

    sed -n '/^begin mrbayes;\$/,/^end;\$/p' ${best_scheme} >> ${bait}.all.nex
    sed -i '/begin mrbayes;/a set autoclose=yes nowarn=yes;' ${bait}.all.nex

    ##################################################
    ### PICK THE FIRST OUTGROUP TAXON PRESENT IN THIS ALIGNMENT
    ##################################################

    OUTGRP=""
    while read -r O; do
        [ -z "\${O}" ] && continue
        if grep -q -w "\${O}" ${bait}.all.nex; then
            OUTGRP="\${O}"
            break
        fi
    done < ${outgroups_list}

    if [ -z "\${OUTGRP}" ]; then
        echo "WARNING: no outgroup taxon found in ${bait}, skipping MrBayes" >&2
        exit 0
    fi

    echo "Using outgroup: \${OUTGRP}"

    sed -i '\$ s/end;//' ${bait}.all.nex

    cat >> ${bait}.all.nex <<EOF
  outgroup \${OUTGRP};
  mcmcp ngen=${params.mb_ngen} samplefreq=${params.mb_samplefreq} mcmcdiagn=yes diagnfreq=${params.mb_diagnfreq} stoprule=yes stopval=${params.mb_stopval};
  mcmc;
  sump;
  sumt conformat=simple;

end;
EOF

    ##################################################
    ### RUN MRBAYES
    ##################################################

    echo "=================================================="
    echo "RUNNING MRBAYES FOR LOCUS ${bait}"
    echo "=================================================="

    set +e
    mb ${bait}.all.nex > ${bait}.log 2>&1
    STATUS=\$?
    set -e

    if [ \${STATUS} -ne 0 ] || [ ! -s "${bait}.all.nex.con.tre" ]; then
        echo "WARNING: MrBayes did not converge / produce a consensus tree for ${bait} (exit=\${STATUS})" >&2
        exit 0
    fi

    ##################################################
    ### CONVERT CONSENSUS TREE TO NEWICK
    ##################################################

    gotree reformat newick --format nexus -i ${bait}.all.nex.con.tre -o ${bait}.all.nex.con.nw
    head -n 1 ${bait}.all.nex.con.nw > ${bait}.BI.treefile

    echo "Finished ${bait} (MrBayes exit=\${STATUS})"
    """
}
