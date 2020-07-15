#!/bin/bash

#$ -V
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

ind=$((${SGE_TASK_ID}-1))


section=0
if [ "$ind" -eq "0" ]
then 
  start_ind=0
else

  num_trees=0
  for i in `seq 0 1 $(($num_windows-1))`
  do
    num_trees=$(( $(cat "relate_chr${c}_${i}.arg" | wc -l) + $num_trees )) #add number of trees in relate_chr${c}_${i}.arg to num_trees
    if [ "$num_trees" -ge "150" ] #if num_trees exceeds 150, increase section by 1
    then
      section=$(($section + 1))
      num_trees=0
    fi

    if [ "$section" -eq "$ind" ] #if section equals $ind, then I found the place to start, which is i+1
    then
      start_ind=$(($i + 1))
      break
    fi
  done

  if [ "$i" -eq $(($num_windows-1)) ] #it 
  then

    exit 0

  fi

fi

if [ $ind -ge $(($num_batched_windows-1)) ] #in this case, I know I need to do this to the last file
then
  end_ind=$num_windows
else
  num_trees=0
  for i in `seq $start_ind 1 $(($num_windows-1))`
  do
    num_trees=$(( $(cat "relate_chr${c}_${i}.arg" | wc -l) + $num_trees ))
    if [ "$num_trees" -ge "150" ] #find $i at which this section ends.
    then
      end_ind=$i
      break
    fi
  done

  if [ "$num_trees" -le "149" ] #if $num_trees is less than 150, it means that $i=$num_windows-1. In this case, end_ind=$num_windows-1
  then
    end_ind=$(($num_windows-1))
  fi

fi

echo $start_ind $end_ind

../bin/./Relate \
    --mode "InferBranchLengths" \
    -m $mu \
    -N $Ne \
    --first_section $start_ind \
    --last_section $end_ind \
    -o ${output} 2>> ${SGE_TASK_ID}_infer_branch_length_c${chunk}.log 

echo "***********************************************"
echo "Finished at: "`date`
echo "***********************************************"
exit 0

