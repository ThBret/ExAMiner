process MAKE_DATABASE_SRA {
    tag "$sample"
    label 'sra_make_db'
    publishDir "${params.sra_db_dir}", mode: 'copy'

    input:
        tuple val(sample), path(r1), path(r2), path(library_conf)

    output:
        tuple val(sample), path("${sample}"), emit: database

    script:
    """
    mkdir -p ${sample}

    mpirun \
        --oversubscribe \
        --mca btl ^vader \
        -np ${task.cpus} \
        SRAssembler_MPI \
        -l ${library_conf} \
        -r . \
        -o ${sample} \
        -P
    """
}
