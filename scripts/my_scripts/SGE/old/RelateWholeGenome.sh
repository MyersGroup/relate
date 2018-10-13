#!/bin/bash

##############
#Execute this in parent directory
#Assuming that all necessary data files are in dir ../

set -u
set -e

if [ $# -le 0 ]
  then
    echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    echo "Not enough arguments supplied. Execute as"
    echo "./RelateWholeGenome [path_to_relate]"
    exit 1;
fi

chromosomes_to_relate=($(cat chromosomes_to_relate.txt))

output="relate"
mv data data_next
mkdir -p data

#rm -rf log/
mkdir -p log

#########################################################

PATH_TO_RELATE=$1

#### copy binaries to directory
mkdir -p bin
cp ${PATH_TO_RELATE}/bin/* bin

#### prepare files
qsub -N convert_from_gp_chr${chromosomes_to_relate[0]} \
     -v c=${chromosomes_to_relate[0]} \
     -wd ${PWD}/data_next/ \
     -o ../log/convert_from_gp.log \
     -e ../log/convert_from_gp.log \
     relate/ConvertFromGP.sh

mu=1.25e-8 #mutation rate
Ne=3e4 #population size
for ((index=0; index < ${#chromosomes_to_relate[@]}; index++))
do

  #move files from data_next to data
  qsub -hold_jid convert_from_gp_chr${chromosomes_to_relate[index]} \
       -N move_files_chr${chromosomes_to_relate[index]} \
       -cwd \
       -o log/move_files.log \
       -e log/move_files.log \
       relate/move_files.sh

  if [ "${chromosomes_to_relate[index]}" -ne "${chromosomes_to_relate[${#chromosomes_to_relate[@]}-1]}" ];
  then
    #### prepare files for next chromosome
    qsub -hold_jid move_files_chr${chromosomes_to_relate[index]} \
         -N convert_from_gp_chr${chromosomes_to_relate[index+1]} \
         -v c=${chromosomes_to_relate[index+1]} \
         -wd ${PWD}/data_next/ \
         -o ../log/convert_from_gp.log \
         -e ../log/convert_from_gp.log \
         relate/ConvertFromGP.sh
  fi

  #### script will wait until chromosome is done, and then submit jobs for the next chromosome
  relate/RelateSGE.sh $PATH_TO_RELATE ${chromosomes_to_relate[index]} $mu $Ne $output
  qsub -cwd -N gzip -v c=${chromosomes_to_relate[index]},output="result/${output}" gzip_results.sh
done

