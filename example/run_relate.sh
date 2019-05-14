#!/bin/bash

#../bin/Relate --mode All \
#  --haps ./data/example.haps.gz \
#  --sample ./data/example.sample.gz \
#  --map ./data/genetic_map.txt \
#  -N 30000 \
#  -m 1.25e-8 \
#  -o example

../scripts/EstimatePopulationSize/EstimatePopulationSize.sh \
  -i example \
  -o example_bypop \
  -m 1.25e-8 \
  --poplabels ./data/example.poplabels \
  --threshold 0 \
  --years_per_gen 28 \
  --num_iter 1 \
  --seed 1

