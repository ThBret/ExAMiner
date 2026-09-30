process PROCESS_READS {
    tag "$sample"
    label 'fastp'
    publishDir "${params.processed_reads_dir}", mode: 'copy', pattern: "*.proc.fastq.gz"
    publishDir "${params.processed_reads_dir}/reports", mode: 'copy', pattern: "*.{html,json}"

    input:
        tuple val(sample), path(r1), path(r2)

    output:
        tuple val(sample), path("${sample}_1.proc.fastq.gz"), path("${sample}_2.proc.fastq.gz"), emit: reads
        tuple path("${sample}.html"), path("${sample}.json"), emit: reports

    script:
    """
    fastp \\
        -i ${r1} -o ${sample}_1.proc.fastq.gz \\
        -I ${r2} -O ${sample}_2.proc.fastq.gz \\
        --dont_overwrite --verbose --detect_adapter_for_pe \\
        --dedup \\
        --dup_calc_accuracy ${params.fastp_dup_acc} \\
        --low_complexity_filter \\
        --complexity_threshold ${params.fastp_complexity} \\
        --trim_poly_g --trim_poly_x \\
        --trim_front1 ${params.fastp_trim_front1} \\
        --trim_front2 ${params.fastp_trim_front2} \\
        --cut_tail \\
        --cut_tail_window_size ${params.fastp_trim_window} \\
        --cut_tail_mean_quality ${params.fastp_trim_quality} \\
        --length_required ${params.fastp_min_length} \\
        --qualified_quality_phred ${params.fastp_qual_phred} \\
        --unqualified_percent_limit ${params.fastp_qual_lim} \\
        --n_base_limit ${params.fastp_base_lim} \\
        --html ${sample}.html --json ${sample}.json \\
        --thread ${task.cpus}
    """
}
