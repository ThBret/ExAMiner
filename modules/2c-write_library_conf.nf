process WRITE_LIBRARY_CONF {
    tag "$sample"
    label 'local_small'
    publishDir "${params.infodir}", mode: 'copy', pattern: "library.*.conf"

    input:
        tuple val(sample), path(r1), path(r2), val(insert_size)

    output:
        tuple val(sample), path(r1), path(r2), path("library.${sample}.conf"), emit: library_conf

    script:
    """
    cat > library.${sample}.conf <<EOF
[LIBRARY]
library_name=${sample}
insert_size=${insert_size}
direction=0
r1=${r1}
r2=${r2}
format=fastq
EOF
    """
}
