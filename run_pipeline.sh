#!/bin/bash
#SBATCH --job-name=beds-test
#SBATCH --output=/scratch/groups/carilee/ngs-pipeline/workflow2/slurm-log.out
#SBATCH --nodes=1
#SBATCH --time=00-12:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=6000
#SBATCH --mail-type=END
#SBATCH --mail-user=habeeb@stanford.edu


set -e
cd $(pwd)
snakemake --cluster-config cluster.json -j 499 \
    --use-singularity \
    --cluster 'sbatch -p {cluster.partition} -t {cluster.time} --mem {cluster.mem} -c {cluster.ncpus} -o {cluster.out}'
