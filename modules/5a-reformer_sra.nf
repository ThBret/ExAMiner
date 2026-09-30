process REFORMER_SRA {
    tag "$bait"
    label 'slurm_small'
    publishDir "${params.reformer_sra_dir}", mode: 'copy'

    input:
        tuple val(bait), val(samples), path(contig_files, stageAs: 'contig_*.fasta'),
              path(bait_stop_fasta, stageAs: 'stop_in.fasta'),
              path(bait_nostop_fasta, stageAs: 'nostop_in.fasta')
    // stageAs to avoid conflicting file names

    output:
        tuple val(bait), path("${bait}"), emit: exonerate_input
        tuple val(bait), path("with.stop/${bait}.fasta"), path("without.stop/${bait}.fasta"), emit: bait_ref

    script:
    def copy_cmds = [samples, contig_files].transpose().collect { s, f ->
        "cp ${f} ${bait}/${s}_all_contigs.fasta"
    }.join('\n    ')
    """
    mkdir -p ${bait} with.stop without.stop
    ${copy_cmds}

    cp ${bait_stop_fasta} with.stop/${bait}.fasta
    cp ${bait_nostop_fasta} without.stop/${bait}.fasta
    """
}
