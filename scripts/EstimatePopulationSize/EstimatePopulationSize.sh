#!/bin/bash

#Treat expanding unset parameters as error 
set -u
#Exit with status 1 if a command throws an error
set -e

if [ $# -le 0 ]
then
  echo ""
  echo "Not enough arguments supplied. Execute as"
  echo "./EstimatePopulationSize.sh"
  echo ""
  echo "--path:      Path to Relate." 
  echo "-i,--input:  Filename of .anc and .mut files without file extension." 
  echo "-o, --output:"
  echo "-m,--mu:     Mutation rate, e.g., 1.25e-8."
  echo "--poplabels: Filename of .poplabels file."
  echo "--threshold: Used to delete trees with less than specified number of mutations. Use 0 to use all trees."
  echo "--num_bins:  "
  echo "--first_chr: "
  echo "--last_chr:  "
  echo ""
  exit 1;
fi

#chr
#num_bins
#threads

######################################################################################################

######################## Read arguments from command line ###########################

POSITIONAL=()
while [[ $# -gt 0 ]]
do
  key="$1"

  case $key in
    --path)
      PATH_TO_RELATE="$2"
      shift # past argument
      shift # past value
      ;;
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

    ${PATH_TO_RELATE}/bin/RelateCoalescentRate \
      --mode ReEstimateBranchLengths \
      --coal ${output}.coal \
      --mrate ${output}_avg.rate \
      --dist ${output}.dist \
      -m 1.25e-8 \
      -i ${output} \
      -o ${output}

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
      -o ${output}_bypop \
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
      -o ${output}_bypop \
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
      ${PATH_TO_RELATE}/bin/RelateCoalescentRate \
        --mode ReEstimateBranchLengths \
        --coal ${output}.coal \
        --mrate ${output}_avg.rate \
        --dist ${output}_chr${chr}.dist \
        -m 1.25e-8 \
        -i ${output}_chr${chr} \
        -o ${output}_chr${chr}
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
      -o ${output}_bypop \
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
      -o ${output}_bypop \
      --first_chr $first_chr \
      --last_chr $last_chr \
      --num_bins ${num_bins} \
      --poplabels ${filename_poplabels}

  fi

fi

#plot results
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
Rscript ${DIR}/plot_population_size.R ${output}


