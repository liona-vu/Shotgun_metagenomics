#!/bin/bash
#Script to download data

#SBATCH --job-name=download_SRR_files
#SBATCH --time=01:00:00
#SBATCH --account=def-cottenie
#SBATCH --cpus-per-task=1
#SBATCH --error=%x_%j_error.txt
#SBATCH --output=%x_%j_output.txt
#SBATCH --mem=4G

# %x is job name, %j is for job id!

#load module
module load sra-toolkit

#make directory and cd into that directory
mkdir data
cd data

#data that will be downloaded and what each sample represents
#SRR8146968 -VOV58 Vegan
#SRR8146970 -Omnivore VOV36
#SRR8146971 - Omnivore VOV28
#SRR8146972 -omnivore VOV26
#SRR8146973 -vegan VOV114
#SRR8146974 -vegan VOV113
#SRR8146975 -omnivore VOV77
#SRR8146976 -omnivore VOV70
#SRR8146977 -vegan VOV29
#SRR8146978 - Vegan VOV20

for i in {68..78};do
  if [[ ${i} = 69 ]]; then #Skipping 69 due to it being vegetarian sample
      continue
  fi
cd SRR81469${i}
fasterq-dump --split-files SRR81469${i}.sra
cd ..
done
