#!/bin/bash
#Script to download data

#load module
module load sra-toolkit

#make directory and cd into that directory
mkdir data
cd data

#data that will be downloaded and what each sample represents
#SRR8146972 -omnivore VOV26
#SRR8146973 -vegan VOV114
#SRR8146974 -vegan VOV113
#SRR8146975 -omnivore VOV77
#SRR8146976 -omnivore VOV70
#SRR8146977 -vegan VOV29

for i in {72..77};
do
echo "Dowloading SRR81469"${i}
prefetch SRR81469${i}
done
