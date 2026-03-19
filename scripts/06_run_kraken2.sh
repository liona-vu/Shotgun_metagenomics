#!/bin/bash

#SBATCH --job-name=kraken2
#SBATCH --time=12:00:00
#SBATCH --account=def-cottenie
#SBATCH --cpus-per-task=32
#SBATCH --error=%x_%j_error.txt
#SBATCH --output=%x_%j_output.txt
#SBATCH --mem=450G

#Use of the kraken2!
module load kraken2/2.1.6

for i in {68..78}; do
  if [[ ${i} = 69 ]]; then
    continue
fi
kraken2 --db database/. \
--confidence 0.15 \ #increase threshold from 0
--threads 32 \
--output SRR81469${i}_trimmed.kraken \
--report SRR81469${i}_trimmed.report \
--paired data/SRR81469${i}/SRR81469${i}_1_trimmed.fastq data/SRR81469${i}/SRR81469${i}_2_trimmed.fastq

done
