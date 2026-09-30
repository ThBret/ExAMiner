process MAPPER {
    tag "$sample"
    label 'mapping'
    publishDir "${params.mapping_out_dir}/6.1-mapping", mode: 'copy', pattern: "*.{bam,bai,txt}"
    publishDir "${params.mapping_out_dir}/refs", mode: 'copy', pattern: "*_loci.fasta*"

    input:
        tuple val(sample), path(bait_fastas), path(r1), path(r2), val(insert_size)

    output:
        tuple val(sample), path("${sample}_loci.fasta"), path("${sample}_loci.fasta.*"), emit: reference
        tuple val(sample), path("${sample}_mapped_mapQ${params.mapping_mapq}_sorted_unique_clipped.bam"),
              path("${sample}_mapped_mapQ${params.mapping_mapq}_sorted_unique_clipped.bam.bai"), emit: bam
        path "${sample}_log*.txt", emit: logs
        path "${sample}_coverage.txt", emit: coverage
        path "${sample}_depth.txt", emit: depth
        path "${sample}_stats.txt", emit: stats
        path "${sample}_flagstats.txt", emit: flagstats

    script:
    """
    # Concatenate all stitched bait loci into one reference
    cat ${bait_fastas} > ${sample}_loci.fasta

    # Index reference for bwa and samtools
    bwa index ${sample}_loci.fasta
    samtools faidx ${sample}_loci.fasta

    # Map reads with bwa mem
    bwa mem -t ${task.cpus} -k ${params.bwa_seed} -U 0 -x intractg \\
        -O ${params.bwa_indel} -E ${params.bwa_gapex} -I ${insert_size} \\
        ${sample}_loci.fasta ${r1} ${r2} -o ${sample}.sam

    # Convert SAM to BAM
    samtools view -@ ${task.cpus} -h -b -t ${sample}_loci.fasta.fai ${sample}.sam -o ${sample}.bam

    # Remove unmapped reads
    samtools view -@ ${task.cpus} -h -b -F 4 ${sample}.bam -o ${sample}_mapped.bam
    samtools view ${sample}_mapped.bam | awk '{print \$3}' | sort | uniq -c > ${sample}_log1_reads_per_locus_mapped.txt

    # MAPQ filtering and sorting
    # Keep only mapped reads higher than a given mapQ and sort the BAM file
    samtools view -@ ${task.cpus} -h -b -q ${params.mapping_mapq} ${sample}.bam | \\
        samtools sort -@ ${task.cpus} -o ${sample}_mapped_mapQ${params.mapping_mapq}_sorted.bam
    samtools view ${sample}_mapped_mapQ${params.mapping_mapq}_sorted.bam | awk '{print \$3}' | sort | uniq -c > ${sample}_log2_reads_per_locus_mapped_mapQ20.txt

    # Check BAM TAGS to remove any reads with multiple hits - this is aligner-specific
    samtools view -@ ${task.cpus} ${sample}_mapped_mapQ${params.mapping_mapq}_sorted.bam | grep -wv "XA:Z" | \\
        samtools view -@ ${task.cpus} -h -b -t ${sample}_loci.fasta.fai - > ${sample}_mapped_mapQ${params.mapping_mapq}_sorted_unique.bam
    samtools view ${sample}_mapped_mapQ${params.mapping_mapq}_sorted_unique.bam | awk '{print \$3}' | sort | uniq -c > ${sample}_log3_reads_per_locus_mapped_mapQ${params.mapping_mapq}_unique.txt

    # Clip overlap of paired-reads to avoid over-estimation of read depth
    bam clipOverlap --in ${sample}_mapped_mapQ${params.mapping_mapq}_sorted_unique.bam \\
        --out ${sample}_mapped_mapQ${params.mapping_mapq}_sorted_unique_clipped.bam
    
    # Index the final BAM file
    samtools index ${sample}_mapped_mapQ${params.mapping_mapq}_sorted_unique_clipped.bam

    # Generate mapping statistics
    samtools coverage ${sample}_mapped_mapQ${params.mapping_mapq}_sorted_unique_clipped.bam > ${sample}_coverage.txt
    samtools depth -@ ${task.cpus} ${sample}_mapped_mapQ${params.mapping_mapq}_sorted_unique_clipped.bam > ${sample}_depth.txt
    samtools stats -@ ${task.cpus} ${sample}_mapped_mapQ${params.mapping_mapq}_sorted_unique_clipped.bam > ${sample}_stats.txt
    samtools flagstats -@ ${task.cpus} ${sample}_mapped_mapQ${params.mapping_mapq}_sorted_unique_clipped.bam > ${sample}_flagstats.txt

    # Remove intermediate files
    rm ${sample}.sam ${sample}.bam ${sample}_mapped_mapQ${params.mapping_mapq}_sorted.bam ${sample}_mapped_mapQ${params.mapping_mapq}_sorted_unique.bam
    """
}
