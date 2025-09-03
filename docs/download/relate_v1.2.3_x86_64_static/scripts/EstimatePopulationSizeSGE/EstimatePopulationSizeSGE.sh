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
  echo "-P:       Assign job to the project."
  echo "-q:       Queue."
  echo "-pe:      Optional. Processor count per node. Default is 'shmem 1'."
  echo "--noanc:           Optional: set --noanc 1 specify to not reinfer anc/mut at the end."
  echo "--pop_of_interest: Optional: Specify the populations for which you want to estimate the population size as a comma separated string. Default is all populations."
  echo "--threshold:       Optional: Used to delete trees with less than specified number of mutations. Use 0 to use all trees. Default: Number of haplotypes."
  echo "--years_per_gen :  Optional: Years per generation. Default is 28."
  echo "--bins:            Optional: Specify epoch bins. Format: lower, upper, stepsize for function c(0,10^seq(lower, upper, stepsize))."
  echo "--num_iter:        Optional: Number of iterations of the algorithm. Default: 5."
  echo "--chr:             Optional: File listing chromosome IDs (one per line). Overrides --first_chr, --last_chr. Assumes input files are indexed by chr, e.g., example_chr1.anc, example_chr1.mut, etc. Specify -i example in this case."
  echo "--first_chr:       Optional: Index of first chr. Assumes that input files are indexed by chr, e.g. example_chr1.anc, example_chr1.mut, etc. Specify -i example in this case. "
  echo "--last_chr:        Optional: Index of last chr."
  echo "--threads:         Optional. Maximum number of threads."
  echo "--seed:            Optional: Random seed for branch lengths estimation"
  echo ""
  exit 1;
fi

#chr
#bins
#threads

PATH_TO_RELATE=$0
PATH_TO_RELATE=$(echo ${PATH_TO_RELATE} | awk -F\scripts/EstimatePopulationSizeSGE/EstimatePopulationSizeSGE.sh '{print $1}')
DIR="${PATH_TO_RELATE}/scripts/EstimatePopulationSize/"

######################################################################################################

######################## Read arguments from command line ###########################

pe="shmem 1"
maxjobs=1
years_per_gen=28
bins="3,7,0.2"

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
    --noanc)
      noanc=1
      shift # past argument
      shift # past value
      ;;
    --bins)
      bins="$2"
      shift # past argument
      shift # past value
      ;;
    --num_iter)
      num_iter="$2"
      shift # past argument
      shift # past value
      ;;
    --chr)
      chr="$2"
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
    -P)
      p="$2"
      shift # past argument
      shift # past value
      ;;
    -q)
      q="$2"
      shift # past argument
      shift # past value
      ;;
    --pe)
      pe="$2"
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
  num_iter=10
else
  if [ $num_iter -le 0 ];
  then
    echo "Error: set num_iter > 0"
    exit 1
  fi
  num_iter_set="FALSE"
fi


echo "********************************"
echo "Parameters passed to script:"
echo "input         = $filename"
echo "poplabels     = $filename_poplabels"
echo "mu            = $mu"
echo "years_per_gen = $years_per_gen"
echo "output        = $output"
echo "num_iter      = $num_iter"
if [ ! -z "${threshold-}" ];
then
  echo "threshold     = $threshold"
fi
if [ ! -z "${pop_of_interest-}" ];
then
  echo "pop_of_interest = $pop_of_interest"
fi
if [ ! -z "${bins-}" ];
then
  echo "bins          = $bins"
fi
if [ ! -z "${first_chr-}" -a ! -z "${last_chr-}" ];
then
  echo "first_chr     = $first_chr"
  echo "last_chr      = $last_chr"
fi
if [ ! -z "${seed-}" ];
then
  echo "seed          = $seed"
fi
echo "Maximum number of threads: $maxjobs" 
echo "********************************"

if [ $filename == $output ];
then
  echo "Please use different names for input and output."
  exit 1;
fi
if [ $maxjobs -le 0 ]
then
  echo "Need positive number of cores."
  exit 1;
fi

######################################################################################################

if [ ! -z "${chr-}" ];
then
  chromosomes=$(cat ${chr}) 
else
  chromosomes=`seq ${first_chr} 1 ${last_chr}`
fi
first_chr=($chromosomes)
filename_chr=${output}.chr
if [ -f ${filename_chr} ]
then
  rm ${filename_chr}
fi
for chr in $chromosomes
do
  echo $chr >> ${filename_chr}
done
tmp=(${chromosomes})
num_chr=${#tmp[@]}

if [ ! -z "${pop_of_interest-}" ];
then

  labels=$(echo ${pop_of_interest} | tr -d ",")

  pop_of_interest=$(echo ${pop_of_interest} | awk '{ gsub(",", ";") ; print $0 }'CC)
  qsub -N extract_subtrees_${output} \
    -sync y \
    -cwd \
    -v PATH_TO_RELATE=${PATH_TO_RELATE},filename=${filename},output=${output},filename_poplabels=${filename_poplabels},pop_of_interest=${pop_of_interest},labels=${labels},filename_chr=${filename_chr} \
    -e ${output}.log \
    -o ${output}.log \
    -t 1-${num_chr} \
    -P $p \
    -q $q \
    -pe shmem 2 \
    ${PATH_TO_RELATE}/scripts/EstimatePopulationSizeSGE/ExtractSubtrees.sh 

  for chr in ${chromosomes}
  do
    cat ${output}_chr${chr}.log  
    rm ${output}_chr${chr}.log  
  done

  mv ${output}_${labels}_chr${first_chr}.poplabels ${output}_${labels}.poplabels
  for chr in `seq $((${first_chr}+1)) 1 ${last_chr}`
  do
    rm ${output}_${labels}_chr${chr}.poplabels 
  done
  filename=${output}_${labels}
  filename_poplabels=${filename}.poplabels

  if [ -z "${threshold-}" ];
  then
    threshold=0.5
  fi

else

  if [ -z "${threshold-}" ];
  then	
    threshold=0.5
  fi

fi

if true
then

qsub -hold_jid extract_subtrees_${output} -N remove_trees_${output} \
  -sync y \
  -cwd \
  -v PATH_TO_RELATE=${PATH_TO_RELATE},filename=${filename},output=${output},threshold=${threshold},filename_chr=${filename_chr} \
  -e ${output}.log \
  -o ${output}.log \
  -t 1-${num_chr} \
  -P $p \
  -q $q \
  -pe ${pe} \
  ${PATH_TO_RELATE}/scripts/EstimatePopulationSizeSGE/RemoveTreesWithFewMutations.sh

for chr in ${chromosomes}
do
  cat ${output}_chr${chr}.log  
  rm ${output}_chr${chr}.log  
done

#remove tmp files of previous runs
foo=$( ls ${output}_chr*_tmp* 2>/dev/null | wc -l)
if [ $foo -ge 1 ]
then
  rm ${output}_chr*_tmp*
fi

bins=$(echo ${bins} | awk '{ gsub(",", ";") ; print $0 }'CC)
qsub -hold_jid remove_trees_${output} -N calc_rates_${output}_0 \
  -cwd \
  -v PATH_TO_RELATE=${PATH_TO_RELATE},filename=${filename},output=${output},years_per_gen=${years_per_gen},bins=${bins} \
  -e ${output}.log \
  -o ${output}.log \
  -P $p \
  -q $q \
  -pe ${pe} \
  ${PATH_TO_RELATE}/scripts/EstimatePopulationSizeSGE/ConstRate.sh 

else
  bins=$(echo ${bins} | awk '{ gsub(",", ";") ; print $0 }'CC)
fi

#repeat iterations of estimating mutation rate, coalescence rates and re-estimating branch lengths
for i in `seq 1 1 ${num_iter}`
do

  qsub -hold_jid calc_rates_${output}_$(($i-1)) -N divide_${output}_${i} \
    -sync y \
    -cwd \
    -v PATH_TO_RELATE=${PATH_TO_RELATE},output=${output},filename_chr=${filename_chr},maxjobs=${maxjobs} \
    -e ${output}.log \
    -o ${output}.log \
    -t 1-${num_chr} \
    -P $p \
    -q $q \
    -pe ${pe} \
    ${PATH_TO_RELATE}/scripts/EstimatePopulationSizeSGE/DivideAncMut.sh 

  for chr in ${chromosomes}
  do

    first_chunk=1
    last_chunk=$(ls ${output}_chr${chr}_tmp_chr*.mut.gz | wc -l)

    if [ -z "${seed-}" ];
    then
      qsub -hold_jid divide_${output}_${i} -N sample_bl_${output}_${i} \
        -cwd \
        -v PATH_TO_RELATE=${PATH_TO_RELATE},output=${output},chr=${chr},mu=${mu} \
        -e ${output}.log \
        -o ${output}.log \
        -t ${first_chunk}-${last_chunk} \
        -P $p \
        -q $q \
        -pe ${pe} \
        ${PATH_TO_RELATE}/scripts/EstimatePopulationSizeSGE/SampleBranchLengths.sh 
    else
      qsub -hold_jid divide_${output}_${i} -N sample_bl_${output}_${i} \
        -cwd \
        -v PATH_TO_RELATE=${PATH_TO_RELATE},output=${output},chr=${chr},mu=${mu},seed=${seed} \
        -e ${output}.log \
        -o ${output}.log \
        -t ${first_chunk}-${last_chunk} \
        -P $p \
        -q $q \
        -pe ${pe} \
        ${PATH_TO_RELATE}/scripts/EstimatePopulationSizeSGE/SampleBranchLengths.sh 
    fi

  done

  qsub -hold_jid sample_bl_${output}_${i} -N combine_${output}_${i} \
    -cwd \
    -v PATH_TO_RELATE=${PATH_TO_RELATE},output=${output},filename_chr=${filename_chr} \
    -e ${output}.log \
    -o ${output}.log \
    -t 1-${num_chr} \
    -P $p \
    -q $q \
    -pe ${pe} \
    ${PATH_TO_RELATE}/scripts/EstimatePopulationSizeSGE/CombineAncMut.sh 

  qsub -hold_jid combine_${output}_${i} -N calc_rates_${output}_${i} \
    -cwd \
    -v PATH_TO_RELATE=${PATH_TO_RELATE},filename=${filename},output=${output},years_per_gen=${years_per_gen},filename_chr=${filename_chr},bins=${bins} \
    -e ${output}.log \
    -o ${output}.log \
    -P $p \
    -q $q \
    -pe ${pe} \
    ${PATH_TO_RELATE}/scripts/EstimatePopulationSizeSGE/CoalRateForTree.sh 	

  cp ${output}.coal ${output}_iter${i}.coal

done

if true
then

#apply pairwise coal
#reestimate branch lengths using coal file
qsub -hold_jid calc_rates_${output}_${num_iter} -N calc_rates_pairwise \
  -cwd \
  -v PATH_TO_RELATE=${PATH_TO_RELATE},filename=${filename},output=${output},years_per_gen=${years_per_gen},filename_chr=${filename_chr},bins=${bins} \
  -t 1-${num_chr} \
  -e ${output}.log \
  -o ${output}.log \
  -P $p \
  -q $q \
  -pe ${pe} \
  ${PATH_TO_RELATE}/scripts/EstimatePopulationSizeSGE/CoalRatePairwise.sh 	

qsub -hold_jid calc_rates_pairwise -N summarize_calc_rates_pairwise \
  -sync y \
  -cwd \
  -v PATH_TO_RELATE=${PATH_TO_RELATE},filename=${filename},output=${output},years_per_gen=${years_per_gen},filename_chr=${filename_chr},bins=${bins} \
  -e ${output}.log \
  -o ${output}.log \
  -P $p \
  -q $q \
  -pe ${pe} \
  ${PATH_TO_RELATE}/scripts/EstimatePopulationSizeSGE/SummarizeCoalRatePairwise.sh 	

for chr in ${chromosomes}
do
  cp ${filename}_chr${chr}.anc.gz ${output}_chr${chr}.anc.gz
  cp ${filename}_chr${chr}.mut.gz ${output}_chr${chr}.mut.gz
done

###########

qsub -hold_jid summarize_calc_rates_pairwise -N divide_${output}_avg \
  -sync y \
  -cwd \
  -v PATH_TO_RELATE=${PATH_TO_RELATE},output=${output},filename_chr=${filename_chr},maxjobs=${maxjobs} \
  -e ${output}.log \
  -o ${output}.log \
  -t 1-${num_chr} \
  -P $p \
  -q $q \
  -pe ${pe} \
  ${PATH_TO_RELATE}/scripts/EstimatePopulationSizeSGE/DivideAncMut.sh 

for chr in ${chromosomes}
do

  first_chunk=1
  last_chunk=$(ls ${output}_chr${chr}_tmp_chr*.mut.gz | wc -l)

  if [ -z "${seed-}" ];
  then
    qsub -hold_jid divide_${output}_${num_iter} -N sample_bl_${output}_avg \
      -cwd \
      -v PATH_TO_RELATE=${PATH_TO_RELATE},output=${output},chr=${chr},mu=${mu} \
      -e ${output}.log \
      -o ${output}.log \
      -t ${first_chunk}-${last_chunk} \
      -P $p \
      -q $q \
      -pe ${pe} \
      ${PATH_TO_RELATE}/scripts/EstimatePopulationSizeSGE/ReEstimateBranchLengths.sh 
  else
    qsub -hold_jid divide_${output}_${num_iter} -N sample_bl_${output}_avg \
      -cwd \
      -v PATH_TO_RELATE=${PATH_TO_RELATE},output=${output},chr=${chr},mu=${mu},seed=${seed} \
      -e ${output}.log \
      -o ${output}.log \
      -t ${first_chunk}-${last_chunk} \
      -P $p \
      -q $q \
      -pe ${pe} \
      ${PATH_TO_RELATE}/scripts/EstimatePopulationSizeSGE/ReEstimateBranchLengths.sh 
  fi

done

qsub -hold_jid sample_bl_${output}_avg -N combine_${output}_avg \
  -cwd \
  -v PATH_TO_RELATE=${PATH_TO_RELATE},output=${output},filename_chr=${filename_chr} \
  -e ${output}.log \
  -o ${output}.log \
  -t 1-${num_chr} \
  -P $p \
  -q $q \
  -pe ${pe} \
  ${PATH_TO_RELATE}/scripts/EstimatePopulationSizeSGE/CombineAncMut.sh 

fi

#plot results
Rscript ${DIR}/plot_population_size.R ${output} ${years_per_gen}


