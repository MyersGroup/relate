#!/bin/bash

#Treat expanding unset parameters as error 
set -u
#Exit with status 1 if a command throws an error
set -e

if [ $# -le 4 ]
  then
    echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    echo "Not enough arguments supplied. Execute as"
    echo "./RelateSGE [path_to_relate] [chr id] [mu] [Ne] [output filename without file extension]"
    exit 1;
fi

######################################################################################################

#### copy binaries to directory
PATH_TO_RELATE=$1
#mkdir -p bin
#cp ${PATH_TO_RELATE}/bin/* bin

#### chromosome index
c=$2

#### create dir for logfiles
#rm -rf log
mkdir -p chr_${c}
mkdir -p chr_${c}/log

output=$5"_chr"$c

#### Set parameters
mu=$3 #mutation rate
Ne=$4 #population size


#### DO NOT CHANGE THE FOLLOWING PARAMETERS
batch_windows=5

######################################################################################################

#### Divide data into chunks of 500k SNPs (the chunks are overlapping to avoid edge effects)
#-sync y causes qsub to wait for the job to complete before exiting. 
qsub -sync y \
     -N make_chunks_chr$c \
     -v c=${c} \
     -wd ${PWD}/chr_${c} \
     -e log/make_chunks.log \
     -o log/make_chunks.log \
     relate/MakeChunks.sh

###########GET NUMBER OF CHUNKS
num_chunks=($(./read_bin.py "chr_${c}/parameters.bin"))
num_chunks=${num_chunks[2]}
prev_chunk=-1

  echo "********************************"
  echo "Number of chunks: "$num_chunks

######################################################################################################
#### Build trees for each chunk

## paint all sequences against each other

for chunk in `seq 0 $(($num_chunks - 1))`;
do

  #GET NUMBER OF WINDOWS
  num_windows=($(./read_bin.py "chr_${c}/parameters_c"${chunk}".bin"))
  num_windows=${num_windows[2]}
  num_batched_windows=$(($num_windows/$batch_windows + 1))

  echo "********************************"
  echo "Parameters passed to scripts:"
  echo "Ne="$Ne "mu="$mu
  echo "Number of windows: "$num_windows

  ## paint all sequences against each other
  #make sure that only 5 paintings exist at a time. Chunk 5 is painted only after chunk 0 is done etc.
  qsub -hold_jid find_equivalent_branches_chr${c}_$(($chunk - 5)) \
       -N paint_chr${c}_${chunk} \
       -wd ${PWD}/chr_${c} \
       -v chunk_index=${chunk} \
       -e log/paint_c${chunk}.log \
       -o log/paint_c${chunk}.log \
       relate/Paint.sh

  ## build tree topologies
  qsub -hold_jid paint_chr${c}_${chunk} \
       -N build_topology_chr${c}_${chunk} \
       -wd ${PWD}/chr_${c} \
       -t 1-$num_batched_windows \
       -v chunk_index=$chunk,output=${output},batch_windows=$batch_windows \
       -e build_c${chunk}.log \
       -o build_c${chunk}.log \
       relate/BuildTopology.sh 

  ## find equivalent branches in adjacent trees 
  qsub -hold_jid build_topology_chr${c}_${chunk} \
       -N find_equivalent_branches_chr${c}_${chunk} \
       -wd ${PWD}/chr_${c}  \
       -v chunk_index=$chunk,output=${output} \
       -e log/find_equivalent_branches_c${chunk}.log \
       -o log/find_equivalent_branches_c${chunk}.log \
       relate/FindEquivalentBranches.sh

  ## infer branch lengths
  qsub -hold_jid find_equivalent_branches_chr${c}_${chunk} \
       -N infer_branch_lengths_chr${c}_${chunk} \
       -wd ${PWD}/chr_${c} \
       -t 1-$num_batched_windows \
       -v Ne=$Ne,mu=$mu,c=${c},chunk_index=$chunk,batch_windows=$batch_windows,output=${output} \
       -e infer_branch_length_c${chunk}.log \
       -o infer_branch_length_c${chunk}.log \
       relate/InferBranchLengths.sh

   ## combine args into one file
  qsub -hold_jid infer_branch_lengths_chr${c}_${chunk} \
       -N combine_args_chr${c} \
       -wd ${PWD}/chr_${c} \
       -v Ne=$Ne,chunk_index=${chunk},c=${c},output=${output} \
       -e log/combine_args_c${chunk}.log \
       -o log/combine_args_c${chunk}.log \
       relate/CombineArgs.sh

  prev_chunk=$chunk

done

############################# Finalize ###########################
#-sync y causes qsub to wait for the job to complete before exiting. 
#finalize results
qsub -sync y \
     -hold_jid combine_args_chr${c} \
     -wd ${PWD}/chr_${c} \
     -N finalize_chr${c} \
     -v output=${output} \
     -e log/combine_args.log \
     -o log/combine_args.log \
     relate/Finalize.sh

#################################################
#clean up directory
for chunk in `seq 0 $(($num_chunks - 1))`;
do
  mv chr_${c}/*.log chr_${c}/log/
done
mv chr_${c}/log chr_${c}/log_chr${c}
tar -cf chr_${c}/log_chr${c}.tar chr_${c}/log_chr${c}
rm -rf chr_${c}/log_chr${c}
mv chr_${c}/log_chr${c}.tar log/
rmdir chr_${c}

qsub -wd ${PWD}/result \
  -N gzip \
  -v c=$c,output="${output}" \
  gzip_results.sh
