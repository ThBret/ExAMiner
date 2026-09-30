process MINING_SRA {
    tag "$sample"
    label 'mining_sra'
    publishDir "${params.sra_mining_dir}", mode: 'copy'

    input:
        tuple val(sample),
              path(database, stageAs: 'database'),
              path(r1),
              path(r2),
              path(baits_list),
              path(baits_stop_fasta),
              path(baits_stop_fai),
              path(baits_stop_id_map),
              path(sra_conf),
              path(library_conf)

    output:
        tuple val(sample), path("${sample}/**"), emit: bait_dirs
        tuple val(sample), path("manifest.tsv"), emit: manifest

    script:
    """
    mkdir -p ${sample} tmp
    : > manifest.tsv

    # Create local SRAssembler library configuration with absolute FASTQ paths
    WORKDIR=\$(pwd)

    sed \\
        -e "s|r1=.*|r1=\${WORKDIR}/${r1.getName()}|" \\
        -e "s|r2=.*|r2=\${WORKDIR}/${r2.getName()}|" \\
        ${library_conf} > library.abs.${sample}.conf

    # Process each bait
    while read -r BAIT; do
        [ -z "\${BAIT}" ] && continue

        OUT=${sample}/\${BAIT}

        # Skip if already assembled
        if [ -s "\${OUT}/all_contigs.fasta" ]; then
            echo "Skipping \${BAIT}: already assembled (\${OUT}/all_contigs.fasta exists)"
            echo -e "\${BAIT}\\t\${OUT}/all_contigs.fasta" >> manifest.tsv
            continue
        fi

        echo "=================================================="
        echo "Starting \${BAIT}"
        date

        RUNDIR=tmp/${sample}.\${BAIT}
        BAIT_WORKDIR=\${RUNDIR}/tmp
        mkdir -p "\${BAIT_WORKDIR}" "\${OUT}"

        # Match bait FASTA by unique ID
        HEADER_ID=\$(awk -F'\t' -v b="\${BAIT}" '\$1==b {print \$2; exit}' ${baits_stop_id_map})

        if [ -z "\${HEADER_ID}" ]; then
            echo "==================================================" >&2
            echo "ERROR: No id-map entry for bait '\${BAIT}'" >&2
            echo "baits_stop_id_map: ${baits_stop_id_map}" >&2
            echo "==================================================" >&2
            continue
        fi

        samtools faidx ${baits_stop_fasta} "\${HEADER_ID}" > "\${RUNDIR}/\${BAIT}.fasta"

        if [ ! -s "\${RUNDIR}/\${BAIT}.fasta" ]; then
            echo "WARNING: faidx extraction empty for '\${BAIT}' (header '\${HEADER_ID}')" >&2
            continue
        fi

        echo "Using bait id: \${BAIT} -> header \${HEADER_ID}"

        cp SRAssembler.conf "\${RUNDIR}/"
        cp library.abs.${sample}.conf "\${RUNDIR}/"

        # Capture absolute paths before changing directory
        WORKDIR="\$(pwd)"
        DB_ABS="\${WORKDIR}/database"
        OUT_ABS="\${WORKDIR}/\${OUT}"
        BAIT_WORKDIR_ABS="\${WORKDIR}/\${BAIT_WORKDIR}"
        QUERY_ABS="\${WORKDIR}/\${RUNDIR}/\${BAIT}.fasta"

        pushd "\${RUNDIR}" > /dev/null

        # ---- signal handling: without this, a SIGTERM from `scancel`
        # ---- (issued when Nextflow is interrupted) gets deferred by bash
        # ---- until the foreground child exits on its own, which can
        # ---- outlast SLURM's kill grace period and leave an orphaned,
        # ---- still-running job behind.

        MPI_PID=""
        cleanup() {
            echo "Caught termination signal - killing srun step (pid \${MPI_PID})" >&2
            [ -n "\${MPI_PID}" ] && kill -TERM "\${MPI_PID}" 2>/dev/null
            sleep 5
            [ -n "\${MPI_PID}" ] && kill -KILL "\${MPI_PID}" 2>/dev/null
            exit 143
        }
        trap cleanup TERM INT

        set +e
        # Run in a subshell so the exported TMPDIR override (needed for a
        # unique PMIx session dir per SLURM job) can never leak into the
        # next bait's iteration of this loop - the while-loop body runs in
        # the same shell process across iterations (no subshell from the
        # `done < file` redirection), so without this scoping an `export`
        # here would remain sticky for every subsequent bait.

        (
            export TMPDIR="\${BAIT_WORKDIR_ABS}/pmix_\${SLURM_JOB_ID:-\$\$}"
            mkdir -p "\${TMPDIR}"
            exec timeout --kill-after=30s ${params.mining_bait_timeout_min}m \\
                mpirun --bind-to none -np ${task.cpus} SRAssembler_MPI </dev/null \\
                -q "\${QUERY_ABS}" \\
                -t protein \\
                -p "\$(basename ${sra_conf})" \\
                -l "\$(basename library.abs.${sample}.conf)" \\
                -r "\${DB_ABS}" \\
                -o "\${OUT_ABS}" \\
                -T "\${BAIT_WORKDIR_ABS}" \\
                -A 1 -k 15:10:45 -S 1 -s drosophila -G 0 \\
                -i 200 -m 200 -M 10000 \\
                -e 0.5 -c 0.8 -n 4 -b 3 -E 3
        ) &

        MPI_PID=\$!
        wait "\${MPI_PID}"
        STATUS=\$?
        trap - TERM INT
        set -e

        # session tmpdir no longer needed once srun exits
        rm -rf "\${BAIT_WORKDIR_ABS}/pmix_\${SLURM_JOB_ID:-\$\$}"
        popd > /dev/null

        echo "Finished \${BAIT} (exit=\${STATUS})"

        if [ \${STATUS} -ne 0 ]; then
            echo "WARNING: SRAssembler failed for '\${BAIT}' (exit=\${STATUS}); continuing with next bait." >&2
            continue
        fi

        if [ -s "\${OUT}/all_contigs.fasta" ]; then
            echo -e "\${BAIT}\\t\${OUT}/all_contigs.fasta" >> manifest.tsv
        else
            echo "WARNING: No all_contigs.fasta produced for '\${BAIT}'." >&2
        fi

    done < ${baits_list}
    """
}
