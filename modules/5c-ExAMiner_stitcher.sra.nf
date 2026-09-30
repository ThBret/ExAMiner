process STITCHER_SRA {
    tag "$bait"
    label 'stitching'
    publishDir "${params.stitcher_out_dir}/${bait}", mode: 'copy', pattern: 'stitcher_out'
    publishDir "${params.stitcher_out_dir}/${bait}", mode: 'copy', pattern: 'mapping_out'
    publishDir "${params.stitcher_out_dir}/${bait}", mode: 'copy', pattern: '*.log'
    publishDir "${params.stitcher_out_dir}/${bait}", mode: 'copy', pattern: 'manifest.tsv'

    input:
        tuple val(bait), path(exonerate_dir), path(bait_nostop_fasta)
        path samples_list
        path outgroups_list

    output:
        tuple val(bait), path("stitcher_out"), emit: stitched
        tuple val(bait), path("mapping_out"), emit: mapping
        // manifest lets main.nf pivot this bait-keyed output into a sample-keyed
        tuple val(bait), path("manifest.tsv"), emit: manifest
        path "${bait}.log", emit: log
        path "mapping_out/tblout/exon_intron_table_${bait}.tsv", optional: true, emit: tblout

    script:
    """
    mkdir -p stitcher_out mapping_out

    # Debugger
    echo "================================================================="
    echo "Running ExAMiner for ${bait}"
    echo "Parameters:"
    echo "  examiner_threshold                   = ${params.examiner_threshold}"
    echo "  examiner_score                       = ${params.examiner_score}"
    echo "  examiner_overlap_mismatch_tolerance  = ${params.examiner_overlap_mismatch_tolerance}"
    echo "  examiner_overlap_conflict_resolution = ${params.examiner_overlap_conflict_resolution}"
    echo "  examiner_missing                     = ${params.examiner_missing}"
    echo "  examiner_buffer                      = ${params.examiner_buffer}"
    echo "  samples_list                         = ${samples_list}"
    echo "  outgroups_list                       = ${outgroups_list}"
    echo "  exonerate_dir                        = ${exonerate_dir}"
    echo "  bait_fasta                           = ${bait_nostop_fasta}"
    echo "  output_dir                           = stitcher_out"
    echo "  mapping_dir                          = mapping_out"
    echo "================================================================="

    echo "Command:"
    echo "python3 ${params.examiner_script} -e ${exonerate_dir} -s ${samples_list} -g ${outgroups_list} -t ${params.examiner_threshold} -S ${params.examiner_score} -u ${params.examiner_missing} -b ${params.examiner_buffer} -o stitcher_out -m mapping_out ${bait_nostop_fasta}"

    # Set error-ignoring state
    set +e
    python3 ${params.examiner_script} \\
        -e ${exonerate_dir} \\
        -s ${samples_list} \\
        -g ${outgroups_list} \\
        -t ${params.examiner_threshold} \\
        -S ${params.examiner_score} \\
        -M ${params.examiner_overlap_mismatch_tolerance} \\
        -C ${params.examiner_overlap_conflict_resolution} \\
        -u ${params.examiner_missing} \\
        -b ${params.examiner_buffer} \\
        -o stitcher_out \\
        -m mapping_out \\
        ${bait_nostop_fasta} > ${bait}.log 2>&1
    STATUS=\$?

    # Unset error-ignoring state
    set -e

    if [ \${STATUS} -ne 0 ]; then
        if grep -q "No candidate samples with reference exons found" "${bait}.log"; then
            echo "WARNING: No candidates for ${bait}, creating empty outputs"
        else
            echo "ERROR: ExAMiner failed unexpectedly for ${bait}" >&2
            cat ${bait}.log >&2
            exit \${STATUS}
        fi
    fi

    : > manifest.tsv
    for d in mapping_out/*/; do
        S=\$(basename "\$d")
        [ "\${S}" = "tblout" ] && continue
        if [ -s "\${d}${bait}.fasta" ]; then
            echo -e "\${S}\\t\${d}${bait}.fasta" >> manifest.tsv
        fi
    done

    echo "Finished ${bait} (ExAMiner exit=\${STATUS})"
    """
}
