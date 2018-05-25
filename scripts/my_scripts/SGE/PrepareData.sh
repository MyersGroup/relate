#!/bin/bash

#$ -V
#$ -j y
#$ -N convert_from_gp
#$ -P myers.prjc -q short.qc

echo "***********************************************"
echo "SGE job ID: "$JOB_ID
echo "SGE task ID: "$SGE_TASK_ID
echo "Run on host: "`hostname`
echo "Operating system: "`uname -s`
echo "Username: "`whoami`
echo "Started at: "`date`
echo "***********************************************"

mkdir chr_${c}

cp ../../1000GP_Phase3/1000GP_Phase3.sample chr_${c}/1000GP_Phase3_chr${c}.sample
cp ../../1000GP_Phase3/genetic_map_chr${c}_combined_b37.txt chr_${c}/
cp ../../human_ancestor_GRCh37_e59/human_ancestor_$c.fa chr_${c}/
cp ../../genome_mask/PilotMask/20140520.chr$c.pilot_mask.fasta chr_${c}/
gunzip -c ../../1000GP_Phase3/1000GP_Phase3_chr$c.hap.gz > chr_${c}/1000GP_Phase3_chr$c.hap
gunzip -c ../../1000GP_Phase3/1000GP_Phase3_chr$c.legend.gz > chr_${c}/1000GP_Phase3_chr$c.legend

pushd chr_${c}

# convert from hap/sample/legend to haps/sample
../../bin/RelateFileFormats \
  --mode ConvertFromHapLegendSample \
  -i 1000GP_Phase3_chr$c \
  --haps chr${c}.haps \
  --sample chr${c}.sample

# flip haps using ancestor
../../bin/RelateFileFormats \
  --mode FlipHapsUsingAncestor \
  --haps chr${c}.haps \
  --sample chr${c}.sample \
  -o chr${c}_ancest \
  --ancestor human_ancestor_${c}.fa 

# Subset individuals
../../bin/RelateFileFormats \
  --mode RemoveSamples \
  --haps chr${c}_ancest.haps \
  --sample chr${c}.sample \
  -i ../remove_ids.txt \
  -o chr${c}_rem

rm chr${c}_ancest.haps
mv chr${c}_rem.sample chr${c}.sample

# filter data using genomic mask
../../bin/RelateFileFormats \
  --mode FilterHapsUsingMask \
  --haps chr${c}_rem.haps \
  --sample chr${c}.sample \
  -o chr${c}_filtered \
  --mask 20140520.chr${c}.pilot_mask.fasta

rm chr${c}_rem.haps
mv chr${c}_filtered.haps chr${c}.haps
mv chr${c}_filtered.dist chr${c}.dist

# generate annotion for SNPs
#../../bin/RelateFileFormats \
#  --mode GenerateSNPAnnotations \
#  --haps chr${c}.haps \
#  --sample chr${c}.sample \
#  --ancestor human_ancestor_${c}.fa \
#  --poplabels 1000GP_Phase3_chr${c}.sample \
#  -o chr${c}

rm 1000GP_Phase3_chr${c}.hap
rm 1000GP_Phase3_chr${c}.legend
rm 1000GP_Phase3_chr${c}.sample
rm human_ancestor_$c.fa
rm 20140520.chr$c.pilot_mask.fasta

gzip chr${c}.haps

popd

echo "***********************************************"
echo "Finished at: "`date`
echo "***********************************************"
exit 0
