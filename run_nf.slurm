#!/bin/bash
#SBATCH --job-name=nextflow_driver
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --partition=long
#SBATCH --output=logs/nextflow_sra_%j.out
#SBATCH --error=logs/nextflow_sra_%j.err

mkdir -p logs

echo "=========================================="
echo "Starting Nextflow pipeline"
echo "Start: $(date)"
echo "Host: $(hostname)"
echo "=========================================="

if [ -d "work" ]; then
    echo "[Nextflow] Found existing work directory, running with -resume"
    nextflow run main.nf -params-file params.yml -resume
else
    echo "[Nextflow] No existing work directory, starting fresh run"
    nextflow run main.nf -params-file params.yml
fi

echo "=========================================="
echo "Finished"
echo "End: $(date)"
echo "=========================================="
