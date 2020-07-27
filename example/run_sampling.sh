#!/bin/bash

../bin/RelateExtract \
	--mode MapMutations \
	--anc example_bypop.anc \
	--mut example_bypop.mut \
	--haps ./data/add.haps \
	--sample ./data/example.sample.gz \
	-o example_mapped

ln -s example_bypop.anc example_mapped.anc

../scripts/SampleBranchLengths/SampleBranchLengths.sh \
	-i example_mapped \
	-o example_sampled \
	--dist example_mapped.dist \
	-m 1.25e-8 \
	--num_samples 5 \
	--first_bp 1000001 \
	--last_bp 1000001 \
	--format b \
	--seed 1 \
	--coal example_bypop.coal

python3 ../scripts/SampleBranchLengths/parse_timeb.py example_sampled.timeb


