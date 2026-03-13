#!/bin/bash

cd data

#to convert SRA file to fastq
for i in {35..40}; 
do 
  cd SRR81469${i}
  echo "Converting SRR81469" ${i}
  fasterq-dump SRR81469${i}.sra
  cd ..
done
