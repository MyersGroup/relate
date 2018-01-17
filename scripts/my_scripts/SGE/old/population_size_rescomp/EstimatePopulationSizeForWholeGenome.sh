#!/bin/bash

#Treat expanding unset parameters as error 
set -u
#Exit with status 1 if a command throws an error
set -e

if [ $# -le 2 ]
  then
    echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    echo "Not enough arguments supplied. Execute as"
    echo "./EstimatePopulationSize [path_to_relate] [pop1] [pop2]"
    exit 1;
fi

PATH_TO_RELATE=$1
pop1=$2
pop2=$3
pop="${pop1}${pop2}"

first_chr=1
last_chr=22

mkdir -p result_${pop}
log="std.log"
#cp -r $PATH_TO_RELATE/bin .

sample=$(ls data/* | grep -e ".sample")
N=$(( 2 * ($(cat $sample | wc -l) - 1) ))

module load R/3.4.0
Rscript get_samples.R $pop1 $pop2
module unload R

############## prepare data ##############

qsub -sync y \
  -hold_jid PrepareFiles \
  -N PrepareFiles \
  -t $first_chr-$last_chr \
  -v pop1=${pop1},pop2=${pop2},pop=${pop} \
  -wd ${PWD}/result_${pop}/ \
  -o \$TASK_ID_std.log \
  -e \$TASK_ID_std.log \
  PrepareFiles.sh

############## Estimate Population Size

## for all chromosomes, do
qsub -sync y \
  -N CoalescentRate_${pop} \
  -t $first_chr-$last_chr \
  -v pop1=${pop1},pop2=${pop2},pop=${pop} \
  -wd ${PWD}/result_${pop}/ \
  -o \$TASK_ID_std.log \
  -e \$TASK_ID_std.log \
  EstimateCoalescentRate.sh

## then do
qsub -sync y \
  -hold_jid CoalescentRate_${pop} \
  -N FinalizeCoalescentRate_${pop} \
  -v pop1=${pop1},pop2=${pop2},pop=${pop},first_chr=$first_chr,last_chr=$last_chr \
  -wd ${PWD}/result_${pop}/ \
  -o std.log \
  -e std.log \
  FinalizeCoalescentRate.sh


############## ReEstimateBranchLengths ################

for i in `seq 1 1 10`
do

  ## for all chromosomes, do
  qsub -sync y \
    -hold_jid FinalizeCoalescentRate_${pop} \
    -N ReEstimateBranchLengths_${pop} \
    -t $first_chr-$last_chr \
    -v pop1=${pop1},pop2=${pop2},pop=${pop} \
    -wd ${PWD}/result_${pop}/ \
    -o \$TASK_ID_std.log \
    -e \$TASK_ID_std.log \
    ReEstimateBranchLengths.sh

  ## for all chromosomes, do
  qsub -sync y \
    -hold_jid ReEstimateBranchLengths_${pop} \
    -N CoalescentRate_${pop} \
    -t $first_chr-$last_chr \
    -v pop1=${pop1},pop2=${pop2},pop=${pop} \
    -wd ${PWD}/result_${pop}/ \
    -o \$TASK_ID_std.log \
    -e \$TASK_ID_std.log \
    EstimateCoalescentRate.sh

  ## then do
  qsub -sync y \
    -hold_jid CoalescentRate_${pop} \
    -N FinalizeCoalescentRate_${pop} \
    -v pop1=${pop1},pop2=${pop2},pop=${pop},first_chr=$first_chr,last_chr=$last_chr \
    -wd ${PWD}/result_${pop}/ \
    -o std.log \
    -e std.log \
    FinalizeCoalescentRate.sh

  cp result_${pop}/population_structure.pdf result_${pop}/population_structure_${i}.pdf
  cp result_${pop}/MutationRateThroughTime_${pop}.pdf result_${pop}/MutationRateThroughTime_${pop}_${i}.pdf

done

if [ "${pop1}" == "${pop2}" ]
then

  ## for all chromosomes, do
  qsub -sync y \
    -hold_jid FinalizeCoalescentRate_${pop} \
    -N move_${pop} \
    -v first_chr=${first_chr},last_chr=${last_chr},pop=${pop} \
    -wd ${PWD}/result_${pop}/ \
    -o \$TASK_ID_std.log \
    -e \$TASK_ID_std.log \
    move.sh

  ## for all chromosomes, do
  qsub -sync y \
    -hold_jid move_${pop} \
    -N ReEstimateBranchLengths_${pop} \
    -t $first_chr-$last_chr \
    -v pop1=${pop1},pop2=${pop2},pop=${pop} \
    -wd ${PWD}/result_${pop}/ \
    -o \$TASK_ID_std.log \
    -e \$TASK_ID_std.log \
    ReEstimateBranchLengths.sh

  ## for all chromosomes, do
  qsub -sync y \
    -hold_jid ReEstimateBranchLengths_${pop} \
    -N CoalescentRate_${pop} \
    -t $first_chr-$last_chr \
    -v pop1=${pop1},pop2=${pop2},pop=${pop} \
    -wd ${PWD}/result_${pop}/ \
    -o \$TASK_ID_std.log \
    -e \$TASK_ID_std.log \
    EstimateCoalescentRate.sh

  ## then do
  qsub -sync y \
    -hold_jid CoalescentRate_${pop} \
    -N FinalizeCoalescentRate_${pop} \
    -v pop1=${pop1},pop2=${pop2},pop=${pop},first_chr=$first_chr,last_chr=$last_chr \
    -wd ${PWD}/result_${pop}/ \
    -o std.log \
    -e std.log \
    FinalizeCoalescentRate.sh

fi

qsub -sync y \
    -hold_jid FinalizeCoalescentRate_${pop} \
    -N gzip_${pop} \
    -wd ${PWD}/result_${pop}/ \
    -o std.log \
    -e std.log \
    gzip.sh
