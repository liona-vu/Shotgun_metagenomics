#!/bin/bash
#Script to check for quality of short read sequencing

#SBATCH --job-name=fastqc
#SBATCH --time=01:30:00
#SBATCH --account=def-cottenie
#SBATCH --cpus-per-task=1
#SBATCH --error=%x_%j_error.txt
#SBATCH --output=%x_%j_output.txt
#SBATCH --mem=10G

#Load module
module load fastqc

cd data

# Perform fastQC to check for quality
for i in {68..78}; do
  if [[ ${i} = 69 ]]; then
      continue
  fi
cd SRR81469${i}
fastqc SRR81469${i}_1.fastq
fastqc SRR81469${i}_2.fastq
cd ..

done
