#!/bin/bash

#Treat expanding unset parameters as error 
set -u
#Exit with status 1 if a command throws an error
set -e

if [ $# -le 0 ]
then
  echo ""
  echo "Not enough arguments supplied. Execute as"
  echo "PATH_TO_RELATE/scripts/EstimatePopulationSize/EstimatePopulationSize.sh"
  echo ""
  echo "-i,--input:   Filename of .anc and .mut files without file extension." 
  echo "-o, --output: Filename of output files without file extension."
  echo "-m,--mu:      Mutation rate, e.g., 1.25e-8."
  echo "--poplabels:  Filename of .poplabels file."
  echo "--pop_of_interest: Optional: Specify the populations for which you want to estimate the population size as a comma separated string. Default is all populations."
  echo "--threshold:       Optional: Used to delete trees with less than specified number of mutations. Use 0 to use all trees. Default: Number of haplotypes."
  echo "--years_per_gen :  Optional: Years per generation. Default is 28."
  echo "--num_bins:        Optional: Number of bins between 1,000ybp and 10,000,000 ybp. Default is 30."
  echo "--num_iter:        Optional: Number of iterations of the algorithm. Default: 5."
  echo "--first_chr:       Optional: Index of first chr. Assumes that input files are indexed by chr, e.g. example_chr1.anc, example_chr1.mut, etc. Specify -i example in this case. "
  echo "--last_chr:        Optional: Index of last chr."
  echo "--coal:            Optional: Initial estimate of coalescent rate, treating all haplotypes as one population. Ignored if num_bins is specified."
  echo "--threads:         Optional. Number of threads used. We currently only parallelize by chromosome."
  echo "--seed:            Optional: Random seed for branch lengths estimation"
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

maxjobs=1
years_per_gen=28

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
    --pop_of_interest)
      pop_of_interest="$2"
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
    --threads)
      maxjobs="$2"
      shift # past argument
      shift # past value
      ;;
    -o|--output)
      output="$2"
      shift # past argument
      shift # past value
      ;;
    --years_per_gen)
      years_per_gen="$2"
      shift # past argument
      shift # past value
      ;;
    --num_bins)
      num_bins="$2"
      shift # past argument
      shift # past value
      ;;
    --num_iter)
      num_iter="$2"
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
    --coal)
      coal="$2"
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

if [ -z "${num_iter-}" ];
then
  num_iter=5
fi

echo "********************************"
echo "Parameters passed to script:"
echo "input         = $filename"
echo "poplabels     = $filename_poplabels"
echo "mu            = $mu"
echo "years_per_gen = $years_per_gen"
echo "output        = $output"
echo "num_iter      = $num_iter"
if [ ! -z "${coal-}" ];
then
  echo "coal          = $coal"
fi
if [ ! -z "${threshold-}" ];
then
  echo "threshold     = $threshold"
fi
if [ ! -z "${pop_of_interest-}" ];
then
  echo "pop_of_interest = $pop_of_interest"
fi
if [ ! -z "${num_bins-}" ];
then
  echo "num_bins      = $num_bins"
fi
if [ ! -z "${first_chr-}" -a ! -z "${last_chr-}" ];
then
  echo "first_chr     = $first_chr"
  echo "last_chr      = $last_chr"
fi
if [ ! -z "${seed-}" ];
then
  echo "seed      = $seed"
fi
echo "********************************"

if [ $filename == $output ];
then
  echo "Please use different names for input and output."
  exit 1;
fi


######################################################################################################

if [ -z "${first_chr-}" -o -z "${last_chr-}" ];
then

  if [ ! -z "${pop_of_interest-}" ];
  then
    
    labels=$(echo ${pop_of_interest} | tr -d ",")
    #delete all trees that have fewer than $threshold mutations
    ${PATH_TO_RELATE}/bin/RelateExtract \
      --mode SubTreesForSubpopulation \
      --poplabels ${filename_poplabels} \
      --pop_of_interest ${pop_of_interest} \
      --anc ${filename}.anc \
      --mut ${filename}.mut \
      -o ${filename}_${labels}

    filename=${filename}_${labels}
    filename_poplabels=${filename}.poplabels

    if [ -z "${threshold-}" ];
    then

      if [ ! -f $filename.anc ]
      then
        threshold=($(zcat $filename.anc.gz | head -1))
      else
        threshold=($(head -1 $filename.anc))
      fi
      threshold=${threshold[1]}

    fi

  else

    if [ -z "${threshold-}" ];
    then

      if [ ! -f $filename.anc ]
      then
        threshold=($(zcat $filename.anc.gz | head -1))
      else
        threshold=($(head -1 $filename.anc))
      fi
      threshold=${threshold[1]}

    fi

  fi


  #delete all trees that have fewer than $threshold mutations
  ${PATH_TO_RELATE}/bin/RelateExtract \
    --mode RemoveTreesWithFewMutations \
    --threshold $threshold \
    --anc ${filename}.anc \
    --mut ${filename}.mut \
    -o ${output} 


  #repeat iterations of estimating mutation rate, coalescence rates and re-estimating branch lengths
  for i in `seq 1 1 ${num_iter}`
  do

    if [ -z "${num_bins-}" ];
    then

      ${PATH_TO_RELATE}/bin/RelateMutationRate \
        --mode Avg \
        --dist ${output}.dist \
        --years_per_gen ${years_per_gen} \
        -i ${output} \
        -o ${output}

      if [ -z "${coal-}" ]
      then
        ${PATH_TO_RELATE}/bin/RelateCoalescentRate \
          --mode EstimatePopulationSize \
          --years_per_gen ${years_per_gen} \
          -i ${output} \
          -o ${output}
      else
        cp $coal ${output}.coal
      fi

    else

      ${PATH_TO_RELATE}/bin/RelateMutationRate \
        --mode Avg \
        --dist ${output}.dist \
        --years_per_gen ${years_per_gen} \
        --num_bins ${num_bins} \
        -i ${output} \
        -o ${output}

      ${PATH_TO_RELATE}/bin/RelateCoalescentRate \
        --mode EstimatePopulationSize \
        --years_per_gen ${years_per_gen} \
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
        -m ${mu} \
        -i ${output} \
        -o ${output}
    else
      ${PATH_TO_RELATE}/bin/RelateCoalescentRate \
        --mode ReEstimateBranchLengths \
        --coal ${output}.coal \
        --mrate ${output}_avg.rate \
        --dist ${output}.dist \
        --seed $seed \
        -m ${mu} \
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
      --years_per_gen ${years_per_gen} \
      -i ${output} \
      -o ${output}

    ${PATH_TO_RELATE}/bin/RelateCoalescentRate \
      --mode EstimatePopulationSize \
      --years_per_gen ${years_per_gen} \
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
      --years_per_gen ${years_per_gen} \
      --num_bins ${num_bins} \
      -i ${output} \
      -o ${output}

    ${PATH_TO_RELATE}/bin/RelateCoalescentRate \
      --mode EstimatePopulationSize \
      --years_per_gen ${years_per_gen} \
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

  if [ ! -z "${pop_of_interest-}" ];
  then
    
    labels=$(echo ${pop_of_interest} | tr -d ",")

    if [ "$maxjobs" -eq 1 ]
    then

      for chr in `seq ${first_chr} 1 ${last_chr}`
      do

        #delete all trees that have fewer than $threshold mutations
        ${PATH_TO_RELATE}/bin/RelateExtract \
          --mode SubTreesForSubpopulation \
          --poplabels ${filename_poplabels} \
          --pop_of_interest ${pop_of_interest} \
          --anc ${filename}_chr${chr}.anc \
          --mut ${filename}_chr${chr}.mut \
          -o ${filename}_${labels}_chr${chr} 

      done

    else

      ExtractSubTrees (){
        chr=$1
        #delete all trees that have fewer than $threshold mutations
        ${PATH_TO_RELATE}/bin/RelateExtract \
          --mode SubTreesForSubpopulation \
          --poplabels ${filename_poplabels} \
          --pop_of_interest ${pop_of_interest} \
          --anc ${filename}_chr${chr}.anc \
          --mut ${filename}_chr${chr}.mut \
          -o ${filename}_${labels}_chr${chr} 2> ${output}_chr${chr}.log  
      }

      parallelize_extract_subtrees () {
        while [ $# -gt 0 ] ; do
          jobcnt=(`jobs -p`)
          if [ ${#jobcnt[@]} -lt $maxjobs ] ; then
            ExtractSubTrees $1 &
            shift
          fi
        done
        wait
      }
      parallelize_extract_subtrees `seq ${first_chr} 1 ${last_chr}`

      for chr in `seq ${first_chr} 1 ${last_chr}`
      do
        cat ${output}_chr${chr}.log  
        rm ${output}_chr${chr}.log  
      done

    fi

    mv ${filename}_${labels}_chr${first_chr}.poplabels ${filename}_${labels}.poplabels
    for chr in `seq $((${first_chr}+1)) 1 ${last_chr}`
    do
      rm ${filename}_${labels}_chr${chr}.poplabels 
    done
    filename=${filename}_${labels}
    filename_poplabels=${filename}.poplabels

    if [ -z "${threshold-}" ];
    then

      if [ ! -f ${filename}_chr${first_chr}.anc ]
      then
        threshold=($(zcat ${filename}_chr${first_chr}.anc.gz | head -1))
      else
        threshold=($(head -1 ${filename}_chr${first_chr}.anc))
      fi
      threshold=${threshold[1]}
 
    fi

  else

    if [ -z "${threshold-}" ];
    then

      if [ ! -f ${filename}_chr${first_chr}.anc ]
      then
        threshold=($(zcat ${filename}_chr${first_chr}.anc.gz | head -1))
      else
        threshold=($(head -1 ${filename}_chr${first_chr}.anc))
      fi
      threshold=${threshold[1]}

    fi

  fi

  if [ "$maxjobs" -eq 1 ]
  then

    for chr in `seq ${first_chr} 1 ${last_chr}`
    do
      #delete all trees that have fewer than $threshold mutations
      ${PATH_TO_RELATE}/bin/RelateExtract \
        --mode RemoveTreesWithFewMutations \
        --threshold $threshold \
        --anc ${filename}_chr${chr}.anc \
        --mut ${filename}_chr${chr}.mut \
        -o ${output}_chr${chr}
    done

  else

    RemoveTreesWithFewMutations (){
      chr=$1
      #delete all trees that have fewer than $threshold mutations
      ${PATH_TO_RELATE}/bin/RelateExtract \
        --mode RemoveTreesWithFewMutations \
        --threshold $threshold \
        --anc ${filename}_chr${chr}.anc \
        --mut ${filename}_chr${chr}.mut \
        -o ${output}_chr${chr} 2> ${output}_chr${chr}.log 
    }

    parallelize_remove_trees () {
      while [ $# -gt 0 ] ; do
        jobcnt=(`jobs -p`)
        if [ ${#jobcnt[@]} -lt $maxjobs ] ; then
          RemoveTreesWithFewMutations $1 &
          shift
        fi
      done
      wait
    }
    parallelize_remove_trees `seq ${first_chr} 1 ${last_chr}`

    for chr in `seq ${first_chr} 1 ${last_chr}`
    do
      cat ${output}_chr${chr}.log  
      rm ${output}_chr${chr}.log  
    done

  fi

  #repeat iterations of estimating mutation rate, coalescence rates and re-estimating branch lengths
  for i in `seq 1 1 ${num_iter}`
  do

    if [ -z "${num_bins-}" ];
    then

      ${PATH_TO_RELATE}/bin/RelateMutationRate \
        --mode Avg \
        --dist ${output} \
        --years_per_gen ${years_per_gen} \
        --first_chr $first_chr \
        --last_chr $last_chr \
        -i ${output} \
        -o ${output}

      if [ -z "${coal-}" ]
      then
        ${PATH_TO_RELATE}/bin/RelateCoalescentRate \
          --mode EstimatePopulationSize \
          --years_per_gen ${years_per_gen} \
          --first_chr $first_chr \
          --last_chr $last_chr \
          -i ${output} \
          -o ${output}
      else
        cp $coal $output.coal
      fi

    else

      ${PATH_TO_RELATE}/bin/RelateMutationRate \
        --mode Avg \
        --dist ${output} \
        --years_per_gen ${years_per_gen} \
        --num_bins ${num_bins} \
        --first_chr $first_chr \
        --last_chr $last_chr \
        -i ${output} \
        -o ${output}

      ${PATH_TO_RELATE}/bin/RelateCoalescentRate \
        --mode EstimatePopulationSize \
        --years_per_gen ${years_per_gen} \
        --num_bins ${num_bins} \
        --first_chr $first_chr \
        --last_chr $last_chr \
        -i ${output} \
        -o ${output}

    fi


    if [ "$maxjobs" -eq 1 ]
    then

      for chr in `seq ${first_chr} 1 ${last_chr}`
      do
        if [ -z "${seed-}" ];
        then
          ${PATH_TO_RELATE}/bin/RelateCoalescentRate \
            --mode ReEstimateBranchLengths \
            --coal ${output}.coal \
            --mrate ${output}_avg.rate \
            --dist ${output}_chr${chr}.dist \
            -m ${mu} \
            -i ${output}_chr${chr} \
            -o ${output}_chr${chr} 
        else
          ${PATH_TO_RELATE}/bin/RelateCoalescentRate \
            --mode ReEstimateBranchLengths \
            --coal ${output}.coal \
            --mrate ${output}_avg.rate \
            --dist ${output}_chr${chr}.dist \
            --seed $seed \
            -m ${mu} \
            -i ${output}_chr${chr} \
            -o ${output}_chr${chr} 
        fi
      done

    else

      ReEstimateBranchLengths (){
     
        chr=$1 
        if [ -z "${seed-}" ];
        then
          ${PATH_TO_RELATE}/bin/RelateCoalescentRate \
            --mode ReEstimateBranchLengths \
            --coal ${output}.coal \
            --mrate ${output}_avg.rate \
            --dist ${output}_chr${chr}.dist \
            -m ${mu} \
            -i ${output}_chr${chr} \
            -o ${output}_chr${chr} 2> ${output}_chr${chr}.log 
        else
          ${PATH_TO_RELATE}/bin/RelateCoalescentRate \
            --mode ReEstimateBranchLengths \
            --coal ${output}.coal \
            --mrate ${output}_avg.rate \
            --dist ${output}_chr${chr}.dist \
            --seed $seed \
            -m ${mu} \
            -i ${output}_chr${chr} \
            -o ${output}_chr${chr} 2> ${output}_chr${chr}.log 
        fi
      
      }

      parallelize_estimating_branchlengths () {
        while [ $# -gt 0 ] ; do
          jobcnt=(`jobs -p`)
          if [ ${#jobcnt[@]} -lt $maxjobs ] ; then
            ReEstimateBranchLengths $1 &
            shift
          fi
        done
        wait
      }
      parallelize_estimating_branchlengths `seq ${first_chr} 1 ${last_chr}`

      for chr in `seq ${first_chr} 1 ${last_chr}`
      do
        cat ${output}_chr${chr}.log  
        rm ${output}_chr${chr}.log  
      done

    fi

  done

  #estimate mutation rate and coalescent rates for a final time

  if [ -z "${num_bins-}" ];
  then

    ${PATH_TO_RELATE}/bin/RelateMutationRate \
      --mode Avg \
      --dist ${output} \
      --years_per_gen ${years_per_gen} \
      --first_chr $first_chr \
      --last_chr $last_chr \
      -i ${output} \
      -o ${output}

    ${PATH_TO_RELATE}/bin/RelateCoalescentRate \
      --mode EstimatePopulationSize \
      --years_per_gen ${years_per_gen} \
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
      --years_per_gen ${years_per_gen} \
      --num_bins ${num_bins} \
      --first_chr $first_chr \
      --last_chr $last_chr \
      -i ${output} \
      -o ${output}

    ${PATH_TO_RELATE}/bin/RelateCoalescentRate \
      --mode EstimatePopulationSize \
      --years_per_gen ${years_per_gen} \
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
      --years_per_gen ${years_per_gen} \
      --num_bins ${num_bins} \
      --poplabels ${filename_poplabels}

  fi

fi

#plot results
DIR="${PATH_TO_RELATE}/scripts/EstimatePopulationSize/"
Rscript ${DIR}/plot_population_size.R ${output} ${years_per_gen}


