#!/bin/bash

##############
#Execute this in dir relateSGE
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

module load python/2.7.10

#chromosomes_to_relate=($(cat chromosomes_to_relate.txt))
chromosomes_to_relate=$(cat chromosomes_to_relate.txt)

output="relate"

#rm -rf log/
mkdir -p log

#########################################################

PATH_TO_RELATE=$1

#### copy binaries to directory
mkdir -p bin
cp ${PATH_TO_RELATE}/bin/* bin

mu=1.25e-8 #mutation rate
Ne=3e4 #population size

maxjobs=2
parallelize () {
  while [ $# -gt 0 ] ; do
    jobcnt=(`jobs -p`)
    if [ ${#jobcnt[@]} -lt $maxjobs ] ; then
      relate/RelateSGE.sh $PATH_TO_RELATE $1 $mu $Ne $output &
      shift
    fi
  done
  wait
}

parallelize ${chromosomes_to_relate}

