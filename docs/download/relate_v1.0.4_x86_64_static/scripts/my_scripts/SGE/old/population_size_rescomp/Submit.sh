#!/bin/bash

populations="ACB ASW BEB CDX CEU CHB CHS CLM ESN FIN GBR GIH GWD IBS ITU JPT KHV LWK MSL MXL PEL PJL PUR STU TSI YRI"
#populations="ACB ASW BEB CDX CEU CHB CHS CLM ESN FIN GIH GWD IBS ITU KHV LWK MSL MXL PEL PJL PUR STU TSI"

for pop in $populations
do
  #echo $pop
  ./EstimatePopulationSizeForWholeGenome.sh ../../relate/ $pop $pop &
  #./ResumePopulationSizeForWholeGenome.sh ../../relate/ $pop $pop &
done

./EstimatePopulationSizeForWholeGenome.sh ../../relate/ GBR YRI &
./EstimatePopulationSizeForWholeGenome.sh ../../relate/ CHB YRI &
./EstimatePopulationSizeForWholeGenome.sh ../../relate/ GBR CHB &
./EstimatePopulationSizeForWholeGenome.sh ../../relate/ GBR FIN &
./EstimatePopulationSizeForWholeGenome.sh ../../relate/ CHB JPT &
./EstimatePopulationSizeForWholeGenome.sh ../../relate/ LWK YRI &

