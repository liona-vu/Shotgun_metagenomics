#!/bin/bash
#script to convert SRA to fastq files

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
