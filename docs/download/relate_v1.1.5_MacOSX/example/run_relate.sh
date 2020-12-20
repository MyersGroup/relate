#!/bin/bash

PATH_TO_RELATE="../"

${PATH_TO_RELATE}/bin/Relate --mode All \
	--haps ./data/example.haps.gz \
	--sample ./data/example.sample.gz \
	--map ./data/genetic_map_GRCh37_chr1.txt \
	-N 30000 \
	-m 1.25e-8 \
	-o example \
	--seed 1

${PATH_TO_RELATE}/scripts/EstimatePopulationSize/EstimatePopulationSize.sh \
	-i example \
	-o example_bypop \
	--noanc 1 \
	-m 1.25e-8 \
	--poplabels ./data/example.poplabels \
	--years_per_gen 28 \
	--seed 1

${PATH_TO_RELATE}/scripts/SampleBranchLengths/ReEstimateBranchLengths.sh \
	-i example \
	-o example_bypop \
	-m 1.25e-8 \
	--coal example_bypop.coal \
	--seed 1

${PATH_TO_RELATE}/scripts/SampleBranchLengths/SampleBranchLengths.sh \
	-i example \
	-o example_bypop_sampled \
	-m 1.25e-8 \
	--coal example_bypop.coal \
	--first_bp 199990000 \
	--last_bp 200010000 \
	--seed 1 \
	--num_samples 100

${PATH_TO_RELATE}/scripts/TreeView/TreeViewSample.sh \
	--haps ./data/example.haps.gz \
	--sample ./data/example.sample.gz \
	--anc example_bypop_sampled.anc \
	--mut example_bypop_sampled.mut \
	--bp_of_interest 200000000 \
	--years_per_gen 28 \
	--poplabels ./data/example.poplabels \
	--dist example_bypop_sampled.dist \
	-o plot
