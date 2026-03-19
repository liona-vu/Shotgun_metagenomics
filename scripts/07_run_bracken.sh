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
for i in {68..78}; do
if [[ ${i} = 69 ]]; then
  continue
fi
  bracken -d database/. \
  -i SRR81469${i}_trimmed.report \
  -o SRR81469${i}.bracken \
  -r 150 \ #short reads were 150bp long
  -l S \ #we want species
  -t 0
done
