#!/bin/bash

#$ -cwd -V
#$ -j y
#$ -P myers.prjc -q short.qc

echo "***********************************************"
echo "SGE job ID: "$JOB_ID
echo "SGE task ID: "$SGE_TASK_ID
echo "Run on host: "`hostname`
echo "Operating system: "`uname -s`
echo "Username: "`whoami`
echo "Started at: "`date`
echo "***********************************************"

chr=${SGE_TASK_ID}

mkdir -p result_chr${chr}
pushd result_chr${chr}

#copy files to working dir
cp ../../../../human_ancestor_GRCh37_e59/human_ancestor_${chr}.fa .
cp ../../../../genome_mask/PilotMask/20140520.chr${chr}.pilot_mask.fasta .
gunzip -c ../../../result/relate_chr${chr}.arg.gz > relate_chr${chr}.arg
gunzip -c ../../../result/relate_chr${chr}.mut.gz > relate_chr${chr}.mut

#create arg and mut files for subpopulation
../../bin/RelateTools --mode CreateArgFileForSubpopulation --pop_of_interest $pop1,$pop2 --samples ../../data/1000GP_Phase3_sub.sample --arg relate_chr${chr}.arg --mut relate_chr${chr}.mut -o relate_chr${chr}_${pop}_pre
#create pos file for subpopulation
../../bin/RelateTools --mode CreatePosFile --mask 20140520.chr${chr}.pilot_mask.fasta --ancestor human_ancestor_${chr}.fa --mut relate_chr${chr}_${pop}_pre.mut -o ../pos_chr${chr}_${pop}.txt 

#remove trees with few mutations
parameters=($(head -1 "relate_chr${chr}_${pop}_pre.arg"))
N=$(( ${parameters[1]} ))
echo $N
../../bin/RelateTools --mode RemoveTreesWithFewMutations --threshold $N --pos ../pos_chr${chr}_${pop}.txt --arg relate_chr${chr}_${pop}_pre.arg --mut relate_chr${chr}_${pop}_pre.mut -o relate_chr${chr}_${pop}_small

#copy back to ..
mv relate_chr${chr}_${pop}_small.arg relate_chr${chr}_${pop}.arg
mv relate_chr${chr}_${pop}_small.mut relate_chr${chr}_${pop}.mut

gzip relate_chr${chr}_${pop}_pre.arg
rm relate_chr${chr}_${pop}_small.pos
rm relate_chr${chr}.arg
rm relate_chr${chr}.mut
rm human_ancestor_${chr}.fa
rm 20140520.chr${chr}.pilot_mask.fasta 

mv * ../

popd

rmdir result_chr${chr}

echo "***********************************************"
echo "Finished at: "`date`
echo "***********************************************"
exit 0

