process DECOMPRESS_READS {
    tag "$sample"
    label 'gunzip'

    input:
        tuple val(sample), path(r1), path(r2)

    output:
        tuple val(sample), path("${sample}_1.proc.fastq"), path("${sample}_2.proc.fastq"), emit: reads

    script:
    """
    gunzip -c ${r1} > ${sample}_1.proc.fastq
    gunzip -c ${r2} > ${sample}_2.proc.fastq
    """
}
