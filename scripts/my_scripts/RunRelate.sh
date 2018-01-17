
Relate \
  --mode "All" \
  --haps ../data/msprime.haps \
  --sample ../data/msprime.sample \
  --map ../data/genetic_map.txt \
  -m 1.25e-8 \
  -N 30000 \
  -o relate
 
mv relate.anc relate_pre.anc
mv relate.mut relate_pre.mut

./EstimatePopulationSize.sh
