/*
 * 9c - GENERATE_DATASETS
 *
 * Final step: build the cross-validation dataset from ML + BI gene trees
 * (datagene.R), copy the whitelisted alignments/trees, root the trees on the
 * outgroups, and derive the metazoa (mzl) BUSCO subset.
 */
process GENERATE_DATASETS {
    tag "crossval_${params.crossval_tag}"
    label 'datasets'
    publishDir "${params.subsets_dir}", mode: 'copy'

    input:
        path bi_trees                                   // collected <bait>.BI.treefile from MRBAYES
        path ml_trees_dir,     stageAs: 'ml_trees'      // published trees dir (<bait>.ML.treefile)
        path alignments_dir,   stageAs: 'alignments_in' // published alignments dir (<bait>.all.fasta)
        path postfilter_stats                           // AMAS stats (ingroup)
        path baits_list
        path outgroups_list
        path busco_table                                // BUSCO full_table.tsv (metazoa_odb10)
        path datagene_script

    output:
        path "crossval_${params.crossval_tag}", emit: crossval

    script:
    def tag     = params.crossval_tag
    def statdir = "crossval_stats_${tag}"
    def out     = "crossval_${tag}"
    """
    echo "=================================================="
    echo "GENERATING CROSSVAL DATASET (${tag}, CID threshold ${params.crossval_cid_thres})"
    echo "=================================================="

    ##################################################
    ### PAIR ML + BI TREES PER LOCUS (DATAGENE INPUT)
    ##################################################

    mkdir -p tmp_inpout

    cp ${postfilter_stats} tmp_inpout/postfilter_stats_ingroup.txt

    N_PAIRED=0
    while read -r B; do
        [ -z "\${B}" ] && continue
        if [ -s "ml_trees/\${B}.ML.treefile" ] && [ -s "\${B}.BI.treefile" ]; then
            cp -L "ml_trees/\${B}.ML.treefile" "\${B}.BI.treefile" tmp_inpout/
            cat "ml_trees/\${B}.ML.treefile" "\${B}.BI.treefile" > "tmp_inpout/\${B}.trees"
            N_PAIRED=\$((N_PAIRED + 1))
        else
            echo "WARNING: ML and/or BI tree missing for \${B}, skipping" >&2
        fi
    done < ${baits_list}

    if [ "\${N_PAIRED}" -eq 0 ]; then
        echo "ERROR: no locus has both an ML and a BI tree" >&2
        exit 1
    fi
    echo "Paired ML+BI trees for \${N_PAIRED} loci"

    ##################################################
    ### GENERATE THE CROSSVAL DATASET
    ##################################################

    Rscript ${datagene_script} ${params.crossval_mode} -r \\
        -i "\$PWD/tmp_inpout/" \\
        -o ${statdir} \\
        -a postfilter_stats_ingroup.txt \\
        -c ${params.crossval_cid_thres}

    WL="tmp_inpout/${statdir}/crossval_whitelist.txt"

    if [ ! -s "\${WL}" ]; then
        echo "ERROR: datagene.R produced no crossval_whitelist.txt" >&2
        exit 1
    fi

    ##################################################
    ### COPY WHITELISTED LOCI
    ##################################################

    mkdir -p ${out}/alignments ${out}/trees
    cp -r tmp_inpout/${statdir} ${out}/

    while read -r W; do
        [ -z "\${W}" ] && continue
        cp -L alignments_in/\${W}.* ${out}/alignments/
        cp -L ml_trees/\${W}.*      ${out}/trees/
        cp -L \${W}.BI.treefile     ${out}/trees/
    done < "\${WL}"

    echo "Whitelisted loci: \$(grep -c . "\${WL}")"

    ##################################################
    ### ROOT THE TREES ON THE OUTGROUPS
    ##################################################

    reroot_tree() {
        if ! gotree reroot outgroup -i "\$1" --strict -l ${outgroups_list} -o "\$2"; then
            echo "WARNING: could not root \$1 (outgroup taxa missing?)" >&2
            rm -f "\$2"
        fi
    }

    while read -r W; do
        [ -z "\${W}" ] && continue
        reroot_tree ${out}/trees/\${W}.ML.treefile ${out}/trees/\${W}.rooted.ML.treefile
        reroot_tree ${out}/trees/\${W}.BI.treefile ${out}/trees/\${W}.rooted.BI.treefile
    done < "\${WL}"

    ##################################################
    ### METAZOA BUSCO SUBSET
    ##################################################

    awk '\$2 == "Complete" {print \$3}' ${busco_table} | sort -u > ${out}/list.busco.mzl.txt

    grep -w -F -f ${out}/list.busco.mzl.txt "\${WL}" \\
        > ${out}/${statdir}/crossval_whitelist_mzlbusco.txt || true

    mkdir -p ${out}/mzlbusco/alignments ${out}/mzlbusco/trees

    while read -r M; do
        [ -z "\${M}" ] && continue
        cp -L ${out}/alignments/\${M}.* ${out}/mzlbusco/alignments/
        cp -L ${out}/trees/\${M}.*      ${out}/mzlbusco/trees/
    done < ${out}/${statdir}/crossval_whitelist_mzlbusco.txt

    echo "BUSCO (mzl) loci in crossval set: \$(grep -c . ${out}/${statdir}/crossval_whitelist_mzlbusco.txt || true)"
    echo "Finished crossval dataset ${tag}"
    """
}
