process CONSENSUS {
    tag "$sample"
    label 'consensus'
    publishDir "${params.consensus_dir}/${params.ucode}/2-stitched.exons", mode: 'copy', pattern: "*.fasta"

    input:
        tuple val(sample), path(bam), path(bai), path(loci_fasta), path(loci_fai)
        path baits_list
        path coordinates

    output:
        tuple val(sample), path("*.fasta"), emit: exons
        tuple val(sample), path("manifest.tsv"), emit: manifest

    script:
    """
    : > manifest.tsv
    mkdir -p exons

    while read -r BAIT; do
        [ -z "\${BAIT}" ] && continue
        grep -q "\${BAIT}" ${loci_fai} || continue

        samtools consensus --ambig -r \${BAIT} -m bayesian -f FASTA -l 0 -a \\
            --min-BQ ${params.samtools_baseq} -d ${params.samtools_depth} \\
            -o \${BAIT}_consensus_d${params.samtools_depth}.fasta ${bam}

        # Get sample and bait-specific exon coordinates
        awk -v S=${sample} -v B=\${BAIT} '(\$1==S && \$2==B && \$3 ~ "exon"){print \$3"\\t"\$4}' \\
            ${coordinates} > \${BAIT}.coordinates.txt

        # Get the exons
        while IFS=\$'\\t' read -r EXONID C; do
            seqkit subseq \${BAIT}_consensus_d${params.samtools_depth}.fasta --region \${C} -w 0 \\
                --update-faidx -o exons/${sample}_\${BAIT}_\${EXONID}.fasta
        done < \${BAIT}.coordinates.txt

        # Stitch the exons together
        N=\$(ls exons/${sample}_\${BAIT}_*.fasta 2>/dev/null | wc -l)
        if [ "\${N}" -eq 1 ]; then
            cp exons/${sample}_\${BAIT}_exon1.fasta \${BAIT}_${sample}.fasta
        elif [ "\${N}" -gt 1 ]; then
            seqkit concat \$(ls exons/${sample}_\${BAIT}_*.fasta | sort -V) -w 0 -o \${BAIT}_${sample}.fasta
        fi

        # Rename the sequences
        if [ -s "\${BAIT}_${sample}.fasta" ]; then
            sed -i '/^>/ s/'\${BAIT}'/'${sample}'/' \${BAIT}_${sample}.fasta
            echo -e "\${BAIT}\\t\${BAIT}_${sample}.fasta" >> manifest.tsv
        fi
    done < ${baits_list}
    """
}
