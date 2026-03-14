#!/bin/bash
#Script to check for quality of short read sequencing

#Load module
module load fastqc

cd data

# Perform fastQC to check for quality
for i in {72..77};
do
  cd SRR81469${i}
  fastqc SRR81469${i}_1.fastq
  fastqc SRR81469${i}_2.fastq
  cd ..
done
