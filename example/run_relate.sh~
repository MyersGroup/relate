#!/bin/bash

if true
then
../bin/Relate --mode All \
  --haps ./data/example.haps.gz \
  --sample ./data/example.sample.gz \
  --map ./data/genetic_map.txt \
  -N 30000 \
  -m 1.25e-8 \
  -o example
fi

if true
then
../scripts/DetectSelection/DetectSelection.sh \
  -i example \
  -o example_bp \
  -m 1.25e-8 \
  --years_per_gen 28 \
  --first_bp 3285308 \
  --last_bp 3285308 

tail example_bp.sele
fi

if true
then

../scripts/EstimatePopulationSize/EstimatePopulationSize.sh \
  -i example \
  -o example_bypop \
  -m 1.25e-8 \
  --poplabels ./data/example.poplabels \
  --threshold 0 \
  --years_per_gen 28 \
  --num_iter 1 \
  --seed 2

  ../scripts/TreeView/TreeView.sh --haps ./data/example.haps.gz --sample ./data/example.sample.gz --poplabels ./data/example.poplabels --anc ./example_bypop.anc --mut ./example_bypop.mut --bp_of_interest 1000000 --years_per_gen 28 -o example

fi

sha1sum example.anc
sha1sum ./before/example.anc

