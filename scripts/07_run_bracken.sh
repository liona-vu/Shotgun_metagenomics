#!/bin/bash

#SBATCH --job-name=bracken
#SBATCH --account=def-cottenie
#SBATCH --cpus-per-task=32
#SBATCH --mem=10G
#SBATCH --time=05:00:00
#SBATCH --error=%x_%j_error.txt
#SBATCH --output=%x_%j_error.txt

module load bracken/3.0
#Loop to run bracken
for i in {72..77};
do
  bracken -d database/. \
  -i SRR81469${i}_trimmed.report \
  -o SRR81469${i}.bracken \
  -r 150 \
  -l S \
  -t 0
done
