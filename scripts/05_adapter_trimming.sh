#!/bin/bash

#SBATCH --job-name=trimming
#SBATCH --time=02:00:00
#SBATCH --account=def-cottenie
#SBATCH --cpus-per-task=1
#SBATCH --error=%x_%j_error.txt
#SBATCH --output=%x_%j_output.txt
#SBATCH --mem=16G

module load fastp/1.0.1

cd data
# trimming for adapters ONLY
for i in {68..78}; do
  if [[ ${i} = 69 ]]; then
    continue
fi
cd SRR81469${i}

fastp --in1 SRR81469${i}_1.fastq \
  --in2 SRR81469${i}_2.fastq \
  --out1 SRR81469${i}_1_trimmed.fastq \
  --out2 SRR81469${i}_2_trimmed.fastq \
  --disable_quality_filtering \
  --disable_length_filtering \
  --detect_adapter_for_pe \
  -h SRR81469${i}.html
  cd ..
done

cd ..
