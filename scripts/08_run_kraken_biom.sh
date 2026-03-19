#!/bin/bash

#to run kraken-biom
kraken-biom \
  SRR8146968_trimmed_bracken_species.report \
  SRR8146970_trimmed_bracken_species.report \
  SRR8146971_trimmed_bracken_species.report \
  SRR8146972_trimmed_bracken_species.report \
  SRR8146973_trimmed_bracken_species.report \
  SRR8146974_trimmed_bracken_species.report \
  SRR8146975_trimmed_bracken_species.report \
  SRR8146976_trimmed_bracken_species.report \
  SRR8146977_trimmed_bracken_species.report \
  SRR8146978_trimmed_bracken_species.report \
  --fmt json \
  -o combined_kraken.biom
