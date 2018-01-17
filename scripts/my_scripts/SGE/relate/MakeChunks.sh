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

#### check if data exists
check_file_existence (){
  if [ ! -f $1 ]; then
      echo "File " + $1 + " not found!"
      exit 1
  fi
}

gunzip -c ../data/chr_${c}/chr${c}.haps.gz > ../data/chr_${c}/chr${c}.haps
check_file_existence ../data/chr_${c}/chr${c}.haps
check_file_existence ../data/chr_${c}/chr${c}.sample
check_file_existence ../data/chr_${c}/chr${c}.dist
check_file_existence ../data/chr_${c}/genetic_map_chr${c}_combined_b37.txt

../bin/./Relate \
  --mode "MakeChunks" \
  --haps ../data/chr_${c}/chr${c}.haps \
  --sample ../data/chr_${c}/chr${c}.sample \
  --map ../data/chr_${c}/genetic_map_chr${c}_combined_b37.txt \
  --dist ../data/chr_${c}/chr${c}.dist \
  --memory 10 2>> log/make_chunks.log

rm ../data/chr_${c}/chr${c}.haps

echo "***********************************************"
echo "Finished at: "`date`
echo "***********************************************"
exit 0
