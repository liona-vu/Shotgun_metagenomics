#!/bin/bash

#download data
mkdir data
cd data

#SRR8146935	- omn
#SRR8146936	- omn
#SRR8146937	- veg
#SRR8146938	- omn
#SRR8146939	- veg
#SRR8146940	- veg

for i in {35..40};
do
echo "Dowloding SRR81469"${i}
prefetch SRR81469${i}
done
