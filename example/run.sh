#!/bin/bash

../scripts/RelateParallel/RelateParallelTopoOnly.sh \
	--haps ./data/example.haps.gz \
	--sample ./data/example.sample.gz \
	--map ./data/genetic_map.txt \
	-N 30000 \
	-m 1.25e-8 \
	-o example2

../bin/Relate --mode All \
	--haps ./data/example.haps.gz \
	--sample ./data/example.sample.gz \
	--map ./data/genetic_map.txt \
	-N 30000 \
	-m 1.25e-8 \
	-o example

../scripts/EstimatePopulationSize/EstimatePopulationSize.sh \
	-i example \
	-o example_bypop \
	-m 1.25e-8 \
	--poplabels ./data/example.poplabels \
	--threshold 0 \
	--years_per_gen 28 \
	--num_iter 1 \
	--seed 1

../scripts/SampleBranchLengths/SampleBranchLengths.sh \
	-i example_bypop \
	-o example_bypop_sampled \
	-m 1.25e-8 \
	--num_samples 3 \
	--first_bp 900000 \
	--last_bp 1000000 \
	--format a \
	--coal example_bypop.coal

../scripts/SampleBranchLengths/SampleBranchLengths.sh \
	-i example2 \
	-o example2_sampled \
	-m 1.25e-8 \
	--num_samples 1 \
	--first_bp 900000 \
	--last_bp 1000000 \
	--format n \
	--coal example_bypop.coal

../scripts/SampleBranchLengths/SampleBranchLengths.sh \
	-i example2 \
	-o example2_sampled \
	-m 1.25e-8 \
	--num_samples 1 \
	--first_bp 900000 \
	--last_bp 1000000 \
	--format b \
	--coal example_bypop.coal
