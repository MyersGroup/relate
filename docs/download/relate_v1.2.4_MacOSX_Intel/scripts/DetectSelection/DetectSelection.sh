#!/bin/bash

#Treat expanding unset parameters as error 
set -u
#Exit with status 1 if a command throws an error
set -e

if [ $# -le 0 ]
then
  echo ""
  echo "Not enough arguments supplied. Execute as"
  echo "PATH_TO_RELATE/scripts/DetectSelection/DetectSelection.sh"
  echo ""
  echo "-i,--input:   Filename of .anc and .mut files without file extension. Needs to contain all SNPs in haps/sample input file." 
  echo "-o, --output: Filename of output files without file extension."
  echo "-m,--mu:      Mutation rate, e.g., 1.25e-8."
  echo "--years_per_gen :  Optional: Years per generation. Default is 28."
  echo "--first_bp:        Optional: First bp position."
  echo "--last_bp:         Optional: Last bp position."
  echo "--coal:            Optional: Estimate of coalescent rate, treating all haplotypes as one population. Not necessary if .anc/.mut files were obtained using EstimatePopulationSize.sh, otherwise strongly recommended."
  echo "--seed:            Optional: Random seed for branch lengths estimation. Only used if --coal is specified."
  echo "--qual:            Optional: Generates ${output}.qual containing heuristics quantifying quality of trees."
  echo ""
  exit 1;
fi

PATH_TO_RELATE=$0
PATH_TO_RELATE=$(echo ${PATH_TO_RELATE} | awk -F\scripts/DetectSelection/DetectSelection.sh '{print $1}')

######################################################################################################

######################## Read arguments from command line ###########################

years_per_gen=28
quality="0"

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
    -m|--mu)
      mu="$2"
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
    --first_bp)
      first_bp="$2"
      shift # past argument
      shift # past value
      ;;
    --last_bp)
      last_bp="$2"
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
    --qual)
      quality="1"
      shift # past argument
      shift # past value
      ;;
    *)    # unknown option
      POSITIONAL+=("$1") # save it in an array for later
      shift # past argument
      ;;
  esac
done


num_iter=0


echo "********************************"
echo "Parameters passed to script:"
echo "input         = $filename"
echo "mu            = $mu"
echo "years_per_gen = $years_per_gen"
echo "output        = $output"

if [ ! -z "${first_bp-}" -a ! -z "${last_bp-}" ];
then
  echo "first_bp      = $first_bp"
  echo "last_bp       = $last_bp"
fi

if [ ! -z "${coal-}" ];
then
  echo "coal          = $coal"
  if [ ! -z "${seed-}" ];
  then
    echo "seed          = $seed"
  fi
else
  if [ ! -z "${seed-}" ];
  then
    echo "seed not used."
  fi
fi
echo "********************************"

#### check if data exists
check_file_existence (){
  if [ ! -f $1 ]; then
    echo "0"
  else
    echo "1"
  fi
}

######################################################################################################

if [ ! -z "${first_bp-}" -a ! -z "${last_bp-}" ];
then

  if [ "${output}" == "${filename}" ]
  then
    echo "Output filename should be different to input filename."
    exit 1
  fi 

  #create output.anc
  ${PATH_TO_RELATE}/bin/RelateExtract \
    --mode AncMutForSubregion \
    --first_bp ${first_bp} \
    --last_bp ${last_bp} \
    --anc ${filename}.anc \
    --mut ${filename}.mut \
    -o ${output}

  num_snps=$(cat ${output}.mut | wc -l)
  if [ ${num_snps} -le 1 ];
  then
    echo "Error: No SNPs included in this region."
    exit 1
  fi

  if [ ! -z "${coal-}" ];
  then

    foo=$(check_file_existence "${filename}.mut")
    num_snps=0
    if [ ${foo} -eq "1" ];
    then
      num_snps=$(cat ${filename}.mut | wc -l)
    else
      foo=$(check_file_existence "${filename}.mut.gz")
      if [ ${foo} -eq "1" ];
      then
        num_snps=$(gunzip -c ${filename}.mut.gz | wc -l)
      fi
    fi

    epochs=$(cat ${coal} | head -2 | tail -1)
    echo -n > ${output}_avg.rate
    for e in ${epochs}
    do
      echo "$e ${mu}" >> ${output}_avg.rate
    done

    if [ ! -z "${seed-}" ];
    then

      #seed set
      ${PATH_TO_RELATE}/bin/RelateCoalescentRate \
        --mode ReEstimateBranchLengths \
        -i ${output} \
        -m ${mu} \
        --dist ${output}.dist \
        --mrate ${output}_avg.rate \
        --coal ${coal} \
        --seed ${seed} \
        -o ${output}

      epochs=$(cat ${output}_avg.rate | awk '{print $1}')
      echo -n > ${output}_avg.rate
      for e in ${epochs}
      do
        echo "$e ${mu}" >> ${output}_avg.rate
      done

      for iter in `seq 1 1 ${num_iter}`
      do
        ${PATH_TO_RELATE}/bin/RelateCoalescentRate \
          --mode ReEstimateBranchLengths \
          -i ${output} \
          -m ${mu} \
          --dist ${output}.dist \
          --mrate ${output}_avg.rate \
          --coal ${coal} \
          --seed ${seed} \
          -o ${output}
      done

    else

      #no seed set
      ${PATH_TO_RELATE}/bin/RelateCoalescentRate \
        --mode ReEstimateBranchLengths \
        -i ${output} \
        -m ${mu} \
        --dist ${output}.dist \
        --mrate ${output}_avg.rate \
        --coal ${coal} \
        -o ${output}
      epochs=$(cat ${output}_avg.rate | awk '{print $1}')
      echo -n > ${output}_avg.rate
      for e in ${epochs}
      do
        echo "$e ${mu}" >> ${output}_avg.rate
      done
      for iter in `seq 1 1 ${num_iter}`
      do
        ${PATH_TO_RELATE}/bin/RelateCoalescentRate \
          --mode ReEstimateBranchLengths \
          -i ${output} \
          -m ${mu} \
          --dist ${output}.dist \
          --mrate ${output}_avg.rate \
          --coal ${coal} \
          -o ${output}
      done

    fi

    rm ${output}.dist
    rm ${output}_avg.rate

    ${PATH_TO_RELATE}/bin/RelateSelection \
      --mode Frequency \
      --years_per_gen ${years_per_gen} \
      -i ${output} \
      -o ${output} 

    ${PATH_TO_RELATE}/bin/RelateSelection \
      --mode Selection \
      -i ${output} \
      -o ${output}

  else

    #only apply selection code
    ${PATH_TO_RELATE}/bin/RelateSelection \
      --mode Frequency \
      --years_per_gen ${years_per_gen} \
      -i ${output} \
      -o ${output} 

    ${PATH_TO_RELATE}/bin/RelateSelection \
      --mode Selection \
      -i ${output} \
      -o ${output}

  fi

else

  if [ ! -z "${coal-}" ];
  then

    echo "!!!!"
    echo "If your genomic region is sufficiently large, we strongly recommend using EstimatePopulationSize.sh first (using --threshold 0)."
    echo "You can then use DetectSelection.sh without specifying the --coal option."
    echo "!!!!"

    if [ ${output} == ${filename} ]
    then
      echo "Output filename should be different to input filename."
      exit 1
    fi 

    foo=$(check_file_existence "${filename}.mut")
    num_snps=0
    if [ ${foo} -eq "1" ];
    then
      num_snps=$(cat ${filename}.mut | wc -l)
    else
      foo=$(check_file_existence "${filename}.mut.gz")
      if [ ${foo} -eq "1" ];
      then
        num_snps=$(gunzip -c ${filename}.mut.gz | wc -l)
      fi
    fi

    epochs=$(cat ${coal} | head -2 | tail -1)
    echo -n > ${output}_avg.rate
    for e in ${epochs}
    do
      echo "$e ${mu}" >> ${output}_avg.rate
    done


    #create output.anc
    cp ${filename}.anc ${output}.anc
    cp ${filename}.mut ${output}.mut

    if [ ! -z "${seed-}" ];
    then

      #seed set
      ${PATH_TO_RELATE}/bin/RelateCoalescentRate \
        --mode ReEstimateBranchLengths \
        -i ${output} \
        -m ${mu} \
        --mrate ${output}_avg.rate \
        --coal ${coal} \
        --seed ${seed} \
        -o ${output}
      
      epochs=$(cat ${output}_avg.rate | awk '{print $1}')
      echo -n > ${output}_avg.rate
      for e in ${epochs}
      do
        echo "$e ${mu}" >> ${output}_avg.rate
      done

      for iter in `seq 1 1 ${num_iter}`
      do
        ${PATH_TO_RELATE}/bin/RelateCoalescentRate \
          --mode ReEstimateBranchLengths \
          -i ${output} \
          -m ${mu} \
          --mrate ${output}_avg.rate \
          --coal ${coal} \
          --seed ${seed} \
          -o ${output}
       done

     else

      #no seed set 
      ${PATH_TO_RELATE}/bin/RelateCoalescentRate \
        --mode ReEstimateBranchLengths \
        -i ${output} \
        -m ${mu} \
        --mrate ${output}_avg.rate \
        --coal ${coal} \
        -o ${output}
      
      epochs=$(cat ${output}_avg.rate | awk '{print $1}')
      echo -n > ${output}_avg.rate
      for e in ${epochs}
      do
        echo "$e ${mu}" >> ${output}_avg.rate
      done
      
      for iter in `seq 1 1 ${num_iter}`
      do
        ${PATH_TO_RELATE}/bin/RelateCoalescentRate \
          --mode ReEstimateBranchLengths \
          -i ${output} \
          -m ${mu} \
          --mrate ${output}_avg.rate \
          --coal ${coal} \
          -o ${output}
      done

    fi

    rm ${output}.dist
    rm ${output}_avg.rate

    ${PATH_TO_RELATE}/bin/RelateSelection \
      --mode Frequency \
      --years_per_gen ${years_per_gen} \
      -i ${output} \
      -o ${output} 

    ${PATH_TO_RELATE}/bin/RelateSelection \
      --mode Selection \
      -i ${output} \
      -o ${output}

  else

    ${PATH_TO_RELATE}/bin/RelateSelection \
      --mode Frequency \
      --years_per_gen ${years_per_gen} \
      -i ${filename} \
      -o ${output} 

    ${PATH_TO_RELATE}/bin/RelateSelection \
      --mode Selection \
      -i ${output} \
      -o ${output}

  fi

fi

#Calculate quality of trees
if [ ${quality} -eq "1" ];
then
  ${PATH_TO_RELATE}/bin/RelateSelection \
    --mode Quality \
    -i ${filename} \
    -o ${output}
fi

