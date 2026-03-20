#!/bin/bash
#script to convert SRA to fastq files

#SBATCH --job-name=convert_SRA
#SBATCH --time=02:00:00
#SBATCH --account=def-cottenie
#SBATCH --cpus-per-task=1
#SBATCH --error=%x_%j_error.txt
#SBATCH --output=%x_%j_output.txt
#SBATCH --mem=10G

#load module
module load sra-toolkit

cd data

#to convert SRA file to fastq
for i in {68..78}; do
  if [[ ${i} = 69 ]]; then
      continue
  fi 
  cd SRR81469${i}
  fasterq-dump --split-files SRR81469${i}.sra
  cd ..
done
