#!/bin/bash

#Treat expanding unset parameters as error 
set -u
#Exit with status 1 if a command throws an error
set -e

if [ $# -le 0 ]
then
  echo ""
  echo "Not enough arguments supplied. Execute as"
  echo "PATH_TO_RELATE/scripts/EstimatePopulationSize.sh"
  echo ""
  echo "-i,--input:  Filename of .anc and .mut files without file extension." 
  echo "-o, --output:"
  echo "-m,--mu:     Mutation rate, e.g., 1.25e-8."
  echo "--poplabels: Filename of .poplabels file."
  echo "--threshold: Optional: Used to delete trees with less than specified number of mutations. Use 0 to use all trees."
  echo "--num_bins:  Optional: Number of bins between 1,000ybp and 10,000,000 ybp. Default is 30."
  echo "--first_chr: Optional: Index of fist chr"
  echo "--last_chr:  Optional: Index of last chr"
  echo "--seed:      Optional: Random seed for branch lengths estimation"
  echo ""
  exit 1;
fi

#chr
#num_bins
#threads

PATH_TO_RELATE=$0
PATH_TO_RELATE=$(echo ${PATH_TO_RELATE} | awk -F\scripts/EstimatePopulationSize/EstimatePopulationSize.sh '{print $1}')

######################################################################################################

######################## Read arguments from command line ###########################

POSITIONAL=()
while [[ $# -gt 0 ]]
do
  key="$1"

  case $key in
    -i|--input)
      filename="$2"
      shift # past argument
      shift # past value
      ;;
    --poplabels)
      filename_poplabels="$2"
      shift # past argument
      shift # past value
      ;;
    -m|--mu)
      mu="$2"
      shift # past argument
      shift # past value
      ;;
    --threshold)
      threshold="$2"
      shift # past argument
      shift # past value
      ;;
    -o|--output)
      output="$2"
      shift # past argument
      shift # past value
      ;;
    --num_bins)
      num_bins="$2"
      shift # past argument
      shift # past value
      ;;
    --first_chr)
      first_chr="$2"
      shift # past argument
      shift # past value
      ;;
    --last_chr)
      last_chr="$2"
      shift # past argument
      shift # past value
      ;;
    --seed)
      seed="$2"
      shift # past argument
      shift # past value
      ;;

    *)    # unknown option
      POSITIONAL+=("$1") # save it in an array for later
      shift # past argument
      ;;
  esac
done

if [ -z "${threshold-}" ];
then
  threshold=($(head -1 example.anc))
  threshold=${threshold[1]}
fi

echo "********************************"
echo "Parameters passed to script:"
echo "input     = $filename"
echo "poplabels = $filename_poplabels"
echo "mu        = $mu"
echo "threshold = $threshold"
echo "output    = $output"
if [ ! -z "${num_bins-}" ];
then
  echo "num_bins  = $num_bins"
fi
if [ ! -z "${first_chr-}" ];
then
  echo "first_chr = $first_chr"
fi
if [ ! -z "${last_chr-}" ];
then
  echo "last_chr  = $last_chr"
fi
if [ ! -z "${seed-}" ];
then
  echo "seed      = $seed"
fi
echo "********************************"


######################################################################################################

if [ -z "${first_chr-}" -o -z "${last_chr-}" ];
then

  #delete all trees that have fewer than $threshold mutations
  ${PATH_TO_RELATE}/bin/RelateExtract \
    --mode RemoveTreesWithFewMutations \
    --threshold $threshold \
    --anc ${filename}.anc \
    --mut ${filename}.mut \
    -o ${output} 


  #repeat iterations of estimating mutation rate, coalescence rates and re-estimating branch lengths
  for i in {1..5}
  do

    if [ -z "${num_bins-}" ];
    then

      ${PATH_TO_RELATE}/bin/RelateMutationRate \
        --mode Avg \
        --dist ${output}.dist \
        -i ${output} \
        -o ${output}

      ${PATH_TO_RELATE}/bin/RelateCoalescentRate \
        --mode EstimatePopulationSize \
        -i ${output} \
        -o ${output}

    else

      ${PATH_TO_RELATE}/bin/RelateMutationRate \
        --mode Avg \
        --dist ${output}.dist \
        --num_bins ${num_bins} \
        -i ${output} \
        -o ${output}

      ${PATH_TO_RELATE}/bin/RelateCoalescentRate \
        --mode EstimatePopulationSize \
        --num_bins ${num_bins} \
        -i ${output} \
        -o ${output}

    fi

    if [ -z "${seed-}" ];
    then
      ${PATH_TO_RELATE}/bin/RelateCoalescentRate \
        --mode ReEstimateBranchLengths \
        --coal ${output}.coal \
        --mrate ${output}_avg.rate \
        --dist ${output}.dist \
        -m 1.25e-8 \
        -i ${output} \
        -o ${output}
    else
      ${PATH_TO_RELATE}/bin/RelateCoalescentRate \
        --mode ReEstimateBranchLengths \
        --coal ${output}.coal \
        --mrate ${output}_avg.rate \
        --dist ${output}.dist \
        --seed $seed \
        -m 1.25e-8 \
        -i ${output} \
        -o ${output}
    fi

  done

  #estimate mutation rate and coalescent rates for a final time

  if [ -z "${num_bins-}" ];
  then

    ${PATH_TO_RELATE}/bin/RelateMutationRate \
      --mode Avg \
      --dist ${output}.dist \
      -i ${output} \
      -o ${output}

    ${PATH_TO_RELATE}/bin/RelateCoalescentRate \
      --mode EstimatePopulationSize \
      -i ${output} \
      -o ${output}

    ${PATH_TO_RELATE}/bin/RelateCoalescentRate \
      --mode FinalizePopulationSize \
      -i ${output} \
      -o ${output} \
      --poplabels ${filename_poplabels}

  else

    ${PATH_TO_RELATE}/bin/RelateMutationRate \
      --mode Avg \
      --dist ${output}.dist \
      --num_bins ${num_bins} \
      -i ${output} \
      -o ${output}

    ${PATH_TO_RELATE}/bin/RelateCoalescentRate \
      --mode EstimatePopulationSize \
      --num_bins ${num_bins} \
      -i ${output} \
      -o ${output}

    ${PATH_TO_RELATE}/bin/RelateCoalescentRate \
      --mode FinalizePopulationSize \
      -i ${output} \
      -o ${output} \
      --num_bins ${num_bins} \
      --poplabels ${filename_poplabels}

  fi

else

  #first_chr and last_chr are specified
  for chr in `seq ${first_chr} 1 ${last_chr}`
  do
    #delete all trees that have fewer than $threshold mutations
    ${PATH_TO_RELATE}/bin/RelateExtract \
      --mode RemoveTreesWithFewMutations \
      --threshold $threshold \
      --anc ${filename}_chr${chr}.anc \
      --mut ${filename}_chr${chr}.mut \
      -o ${output} 
  done

  #repeat iterations of estimating mutation rate, coalescence rates and re-estimating branch lengths
  for i in {1..5}
  do

    if [ -z "${num_bins-}" ];
    then

      ${PATH_TO_RELATE}/bin/RelateMutationRate \
        --mode Avg \
        --dist ${output} \
        --first_chr $first_chr \
        --last_chr $last_chr \
        -i ${output} \
        -o ${output}

      ${PATH_TO_RELATE}/bin/RelateCoalescentRate \
        --mode EstimatePopulationSize \
        --first_chr $first_chr \
        --last_chr $last_chr \
        -i ${output} \
        -o ${output}

    else

      ${PATH_TO_RELATE}/bin/RelateMutationRate \
        --mode Avg \
        --dist ${output} \
        --num_bins ${num_bins} \
        --first_chr $first_chr \
        --last_chr $last_chr \
        -i ${output} \
        -o ${output}

      ${PATH_TO_RELATE}/bin/RelateCoalescentRate \
        --mode EstimatePopulationSize \
        --num_bins ${num_bins} \
        --first_chr $first_chr \
        --last_chr $last_chr \
        -i ${output} \
        -o ${output}

    fi

    for chr in `seq ${first_chr} 1 ${last_chr}`
    do
      if [ -z "${seed-}" ];
      then
        ${PATH_TO_RELATE}/bin/RelateCoalescentRate \
          --mode ReEstimateBranchLengths \
          --coal ${output}.coal \
          --mrate ${output}_avg.rate \
          --dist ${output}.dist \
          -m 1.25e-8 \
          -i ${output} \
          -o ${output}
      else
        ${PATH_TO_RELATE}/bin/RelateCoalescentRate \
          --mode ReEstimateBranchLengths \
          --coal ${output}.coal \
          --mrate ${output}_avg.rate \
          --dist ${output}.dist \
          --seed $seed \
          -m 1.25e-8 \
          -i ${output} \
          -o ${output}
      fi
    done

  done

  #estimate mutation rate and coalescent rates for a final time

  if [ -z "${num_bins-}" ];
  then

    ${PATH_TO_RELATE}/bin/RelateMutationRate \
      --mode Avg \
      --dist ${output} \
      --first_chr $first_chr \
      --last_chr $last_chr \
      -i ${output} \
      -o ${output}

    ${PATH_TO_RELATE}/bin/RelateCoalescentRate \
      --mode EstimatePopulationSize \
      --first_chr $first_chr \
      --last_chr $last_chr \
      -i ${output} \
      -o ${output}

    ${PATH_TO_RELATE}/bin/RelateCoalescentRate \
      --mode FinalizePopulationSize \
      --first_chr $first_chr \
      --last_chr $last_chr \
      -i ${output} \
      -o ${output} \
      --poplabels ${filename_poplabels}

  else

    ${PATH_TO_RELATE}/bin/RelateMutationRate \
      --mode Avg \
      --dist ${output} \
      --num_bins ${num_bins} \
      --first_chr $first_chr \
      --last_chr $last_chr \
      -i ${output} \
      -o ${output}

    ${PATH_TO_RELATE}/bin/RelateCoalescentRate \
      --mode EstimatePopulationSize \
      --num_bins ${num_bins} \
      --first_chr $first_chr \
      --last_chr $last_chr \
      -i ${output} \
      -o ${output}

    ${PATH_TO_RELATE}/bin/RelateCoalescentRate \
      --mode FinalizePopulationSize \
      -i ${output} \
      -o ${output} \
      --first_chr $first_chr \
      --last_chr $last_chr \
      --num_bins ${num_bins} \
      --poplabels ${filename_poplabels}

  fi

fi

#plot results
DIR="${PATH_TO_RELATE}/scripts/EstimatePopulationSize/"
Rscript ${DIR}/plot_population_size.R ${output}


