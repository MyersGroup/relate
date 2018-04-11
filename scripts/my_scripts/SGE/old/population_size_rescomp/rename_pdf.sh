#!/bin/bash

mkdir plots

populations="ACB ASW BEB CDX CEU CHB CHS CLM ESN FIN GBR GIH GWD IBS ITU JPT KHV LWK MSL MXL PEL PJL PUR STU TSI YRI"

for p in $population
do
  cp result_${p}${p}/population_structure.pdf.gz plots/population_structure_${p}${p}.pdf.gz
done

cp result_GBRYRI/population_structure.pdf.gz plots/population_structure_GBRYRI.pdf.gz
cp result_CHBYRI/population_structure.pdf.gz plots/population_structure_CHBYRI.pdf.gz
cp result_GBRCHB/population_structure.pdf.gz plots/population_structure_GBRCHB.pdf.gz
cp result_GBRFIN/population_structure.pdf.gz plots/population_structure_GBRFIN.pdf.gz
cp result_CHBJPT/population_structure.pdf.gz plots/population_structure_CHBJPT.pdf.gz
cp result_LWKYRI/population_structure.pdf.gz plots/population_structure_LWKYRI.pdf.gz
