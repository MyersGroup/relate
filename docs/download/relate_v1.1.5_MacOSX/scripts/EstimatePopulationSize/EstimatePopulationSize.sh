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
PATH_TO_RELATE=$(echo ${PATH_TO_RELATE} | awk -F\scripts/EstimatePopulationSize/EstimatePopulationSize.sh '{print $1}')
DIR="${PATH_TO_RELATE}/scripts/EstimatePopulationSize/"

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

# if only one core available
if [ "$maxjobs" -eq 1 ]
then


  if [ -z "${first_chr-}" -o -z "${last_chr-}" ] && [ -z "${chr-}" ];
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
        -o ${output}_${labels}

      filename=${output}_${labels}
      filename_poplabels=${filename}.poplabels

      if [ -z "${threshold-}" ];
      then

        #if [ ! -f $filename.anc ]
        #then
        #  threshold=($(gunzip -c $filename.anc.gz | head -1))
        #else
        #  threshold=($(head -1 $filename.anc))
        #fi
        #threshold=${threshold[1]}
        threshold=0.5

      fi

    else

      if [ -z "${threshold-}" ];
      then

        #if [ ! -f $filename.anc ]
        #then
        #  threshold=($(gunzip -c $filename.anc.gz | head -1))
        #else
        #  threshold=($(head -1 $filename.anc))
        #fi
        #threshold=${threshold[1]}
        threshold=0.5

      fi

    fi

    #delete all trees that have fewer than $threshold mutations
    ${PATH_TO_RELATE}/bin/RelateExtract \
      --mode RemoveTreesWithFewMutations \
      --threshold $threshold \
      --anc ${filename}.anc \
      --mut ${filename}.mut \
      -o ${output} 

   if [ -z "${bins-}" ];
    then

      ${PATH_TO_RELATE}/bin/RelateCoalescentRate \
        --mode CoalRateForTree \
        --years_per_gen ${years_per_gen} \
        --dist ${output}.dist \
        -i ${output} \
        -o ${output}

    else

      ${PATH_TO_RELATE}/bin/RelateCoalescentRate \
        --mode CoalRateForTree \
        --years_per_gen ${years_per_gen} \
        --bins ${bins} \
        --dist ${output}.dist \
        -i ${output} \
        -o ${output}

    fi

    #repeat iterations of estimating mutation rate, coalescence rates and re-estimating branch lengths
    for i in `seq 1 1 ${num_iter}`
    do

      if [ -z "${seed-}" ];
      then

        ${PATH_TO_RELATE}/bin/RelateCoalescentRate \
          --mode SampleBranchLengths \
          --coal ${output}.coal \
          --dist ${output}.dist \
          --num_samples 1 \
          -m ${mu} \
          -i ${output} \
          -o ${output}

      else

        ${PATH_TO_RELATE}/bin/RelateCoalescentRate \
          --mode SampleBranchLengths \
          --coal ${output}.coal \
          --dist ${output}.dist \
          --num_samples 1 \
					--seed $(($seed + $i)) \
          -m ${mu} \
          -i ${output} \
          -o ${output}

      fi

      if [ -z "${bins-}" ];
      then

        ${PATH_TO_RELATE}/bin/RelateCoalescentRate \
          --mode CoalRateForTree \
          --years_per_gen ${years_per_gen} \
          --coal ${output}.coal \
          --dist ${output}.dist \
          -i ${output} \
          -o ${output}

      else

        ${PATH_TO_RELATE}/bin/RelateCoalescentRate \
          --mode CoalRateForTree \
          --years_per_gen ${years_per_gen} \
          --bins ${bins} \
          --coal ${output}.coal \
          --dist ${output}.dist \
          -i ${output} \
          -o ${output}

      fi 

    done

    if [ -z "${bins-}" ];
    then
      ${PATH_TO_RELATE}/bin/RelateCoalescentRate \
        --mode EstimatePopulationSize \
        --years_per_gen ${years_per_gen} \
        --dist ${output}.dist \
        -i ${output} \
        -o ${output}.pairwise
    else
      ${PATH_TO_RELATE}/bin/RelateCoalescentRate \
        --mode EstimatePopulationSize \
        --years_per_gen ${years_per_gen} \
        --bins ${bins} \
        --dist ${output}.dist \
        -i ${output} \
        -o ${output}.pairwise
    fi

    ${PATH_TO_RELATE}/bin/RelateMutationRate \
      --mode Avg \
      --dist ${output}.dist \
      --years_per_gen ${years_per_gen} \
      -i ${output} \
      -o ${output}

    #estimate mutation rate and coalescent rates for a final time
		if [ -z ${noanc-} ];
		then
			if [ -z "${seed-}" ];
			then

				${PATH_TO_RELATE}/bin/RelateCoalescentRate \
					--mode ReEstimateBranchLengths \
					--coal ${output}.coal \
					--dist ${output}.dist \
					-m ${mu} \
					-i ${filename} \
					-o ${output}

			else

				${PATH_TO_RELATE}/bin/RelateCoalescentRate \
					--mode ReEstimateBranchLengths \
					--coal ${output}.coal \
					--dist ${output}.dist \
					--seed $(($seed + ${num_iter} + 1)) \
					-m ${mu} \
					-i ${filename} \
					-o ${output}

			fi
	  else
			rm ${output}.anc
			rm ${output}.mut
			rm ${output}.dist
		fi

  else

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

    if [ ! -z "${pop_of_interest-}" ];
    then

      labels=$(echo ${pop_of_interest} | tr -d ",")

      for chr in $chromosomes
      do

        #delete all trees that have fewer than $threshold mutations
        ${PATH_TO_RELATE}/bin/RelateExtract \
          --mode SubTreesForSubpopulation \
          --poplabels ${filename_poplabels} \
          --pop_of_interest ${pop_of_interest} \
          --anc ${filename}_chr${chr}.anc \
          --mut ${filename}_chr${chr}.mut \
          -o ${output}_${labels}_chr${chr} 

        done

        mv ${output}_${labels}_chr${first_chr}.poplabels ${output}_${labels}.poplabels
        for chr in $chromosomes
        do
          rm ${output}_${labels}_chr${chr}.poplabels 
        done
        filename=${output}_${labels}
        filename_poplabels=${filename}.poplabels

        if [ -z "${threshold-}" ];
        then

          #if [ ! -f ${filename}_chr${first_chr}.anc ]
          #then
          #  threshold=($(gunzip -c ${filename}_chr${first_chr}.anc.gz | head -1))
          #else
          #  threshold=($(head -1 ${filename}_chr${first_chr}.anc))
          #fi
          #threshold=${threshold[1]}
          threshold=0.5

        fi

      else

        if [ -z "${threshold-}" ];
        then

          #if [ ! -f ${filename}_chr${first_chr}.anc ]
          #then
          #  threshold=($(gunzip -c ${filename}_chr${first_chr}.anc.gz | head -1))
          #else
          #  threshold=($(head -1 ${filename}_chr${first_chr}.anc))
          #fi
          #threshold=${threshold[1]}
          threshold=0.5

        fi

    fi

    for chr in ${chromosomes}
    do

      #delete all trees that have fewer than $threshold mutations
      ${PATH_TO_RELATE}/bin/RelateExtract \
        --mode RemoveTreesWithFewMutations \
        --threshold $threshold \
        --anc ${filename}_chr${chr}.anc \
        --mut ${filename}_chr${chr}.mut \
        -o ${output}_chr${chr}

    done

    if [ -z "${bins-}" ];
    then

      ${PATH_TO_RELATE}/bin/RelateCoalescentRate \
        --mode CoalRateForTree \
        --years_per_gen ${years_per_gen} \
        --dist ${output} \
        --chr ${filename_chr} \
        -i ${output} \
        -o ${output}

    else

      ${PATH_TO_RELATE}/bin/RelateCoalescentRate \
        --mode CoalRateForTree \
        --years_per_gen ${years_per_gen} \
        --bins ${bins} \
        --chr ${filename_chr} \
        --dist ${output} \
        -i ${output} \
        -o ${output}

    fi

    #repeat iterations of estimating mutation rate, coalescence rates and re-estimating branch lengths
    for i in `seq 1 1 ${num_iter}`
    do

      for chr in ${chromosomes}
      do
        if [ -z "${seed-}" ];
        then

          ${PATH_TO_RELATE}/bin/RelateCoalescentRate \
            --mode SampleBranchLengths \
            --coal ${output}.coal \
            --dist ${output}_chr${chr}.dist \
            --num_samples 1 \
            -m ${mu} \
            -i ${output}_chr${chr} \
            -o ${output}_chr${chr}

        else

          ${PATH_TO_RELATE}/bin/RelateCoalescentRate \
            --mode SampleBranchLengths \
            --coal ${output}.coal \
            --dist ${output}_chr${chr}.dist \
            --num_samples 1 \
						--seed $(( $seed + $i )) \
            -m ${mu} \
            -i ${output}_chr${chr} \
            -o ${output}_chr${chr} 

        fi
      done

      if [ -z "${bins-}" ];
      then

        ${PATH_TO_RELATE}/bin/RelateCoalescentRate \
          --mode CoalRateForTree \
          --years_per_gen ${years_per_gen} \
          --chr ${filename_chr} \
          --coal ${output}.coal \
          --dist ${output} \
          -i ${output} \
          -o ${output} 

      else

        ${PATH_TO_RELATE}/bin/RelateCoalescentRate \
          --mode CoalRateForTree \
          --years_per_gen ${years_per_gen} \
          --bins ${bins} \
          --chr ${filename_chr} \
          --coal ${output}.coal \
          --dist ${output} \
          -i ${output} \
          -o ${output}

      fi

    done

    if [ -z "${bins-}" ];
    then
      ${PATH_TO_RELATE}/bin/RelateCoalescentRate \
        --mode EstimatePopulationSize \
        --years_per_gen ${years_per_gen} \
        --chr ${filename_chr} \
        --dist ${output} \
        -i ${output} \
        -o ${output}.pairwise
    else
      ${PATH_TO_RELATE}/bin/RelateCoalescentRate \
        --mode EstimatePopulationSize \
        --years_per_gen ${years_per_gen} \
        --chr ${filename_chr} \
        --bins ${bins} \
        --dist ${output} \
        -i ${output} \
        -o ${output}.pairwise
    fi

    ${PATH_TO_RELATE}/bin/RelateMutationRate \
      --mode Avg \
      --dist ${output} \
      --years_per_gen ${years_per_gen} \
      --chr ${filename_chr} \
      -i ${output} \
      -o ${output}

    #estimate mutation rate and coalescent rates for a final time
		if [ -z ${noanc-} ]
		then
			for chr in ${chromosomes}
			do
				if [ -z "${seed-}" ];
				then

					${PATH_TO_RELATE}/bin/RelateCoalescentRate \
						--mode ReEstimateBranchLengths \
						--coal ${output}.coal \
						--dist ${output}_chr${chr}.dist \
						-m ${mu} \
						-i ${filename}_chr${chr} \
						-o ${output}_chr${chr} 

				else

					${PATH_TO_RELATE}/bin/RelateCoalescentRate \
						--mode ReEstimateBranchLengths \
						--coal ${output}.coal \
						--dist ${output}_chr${chr}.dist \
						--seed $(( $seed + ${num_iter}+1 )) \
						-m ${mu} \
						-i ${filename}_chr${chr} \
						-o ${output}_chr${chr} 

        fi
      done
  	else
			for chr in ${chromosomes}
			do
				rm ${output}_chr${chr}.anc
				rm ${output}_chr${chr}.mut
				rm ${output}_chr${chr}.dist
			done
		fi

  fi

else

  if [ -z "${first_chr-}" -o -z "${last_chr-}" ] && [ -z "${chr-}" ];
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
        -o ${output}_${labels}

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


    #delete all trees that have fewer than $threshold mutations
    ${PATH_TO_RELATE}/bin/RelateExtract \
      --mode RemoveTreesWithFewMutations \
      --threshold $threshold \
      --anc ${filename}.anc \
      --mut ${filename}.mut \
      -o ${output} 

    gzip ${output}.anc
    gzip ${output}.mut

    #remove tmp files of previous runs
    foo=$(ls ${output}_tmp* 2>/dev/null | wc -l)
    if [ $foo -ge 1 ];
    then
      rm ${output}_tmp*
    fi

    if [ -z "${bins-}" ];
    then

      ${PATH_TO_RELATE}/bin/RelateCoalescentRate \
        --mode CoalRateForTree \
        --years_per_gen ${years_per_gen} \
        --dist ${output}.dist \
        -i ${output} \
        -o ${output}

    else

      ${PATH_TO_RELATE}/bin/RelateCoalescentRate \
        --mode CoalRateForTree \
        --years_per_gen ${years_per_gen} \
        --bins ${bins} \
        --dist ${output}.dist \
        -i ${output} \
        -o ${output}

    fi

    #repeat iterations of estimating mutation rate, coalescence rates and re-estimating branch lengths
    for i in `seq 1 1 ${num_iter}`
    do

      ${PATH_TO_RELATE}/bin/RelateExtract \
        --mode DivideAncMut \
        --threads $maxjobs \
        --anc ${output}.anc \
        --mut ${output}.mut \
        -o ${output}_tmp

      first_chunk=0
      last_chunk=$(ls ${output}_tmp_chr*.mut.gz | wc -l)
      last_chunk=$((${last_chunk} - 1)) 

      ReEstimateBranchLengths (){

        chunk=$1 
        if [ -z "${seed-}" ];
        then

          ${PATH_TO_RELATE}/bin/RelateCoalescentRate \
            --mode SampleBranchLengths \
            --coal ${output}.coal \
            --dist ${output}.dist \
            --num_samples 1 \
            -m ${mu} \
            -i ${output}_tmp_chr${chunk} \
            -o ${output}_tmp_chr${chunk} 2> ${output}_tmp_chr${chunk}.log 

        else

          ${PATH_TO_RELATE}/bin/RelateCoalescentRate \
            --mode SampleBranchLengths \
            --coal ${output}.coal \
            --dist ${output}.dist \
            --num_samples 1 \
						--seed $(( $seed + $i )) \
            -m ${mu} \
            -i ${output}_tmp_chr${chunk} \
            -o ${output}_tmp_chr${chunk} 2> ${output}_tmp_chr${chunk}.log 

        fi

        rm ${output}_tmp_chr${chunk}.anc.gz
        rm ${output}_tmp_chr${chunk}.mut.gz

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
      parallelize_estimating_branchlengths `seq ${first_chunk} 1 ${last_chunk}`

			for chr in `seq ${first_chunk} 1 ${last_chunk}`
			do
				cat ${output}_tmp_chr${chr}.log  
				rm ${output}_tmp_chr${chr}.log  
			done

			${PATH_TO_RELATE}/bin/RelateExtract \
				--mode CombineAncMut \
				-o ${output}_tmp

			mv ${output}_tmp.anc.gz ${output}.anc.gz
			mv ${output}_tmp.mut.gz ${output}.mut.gz

      if [ -z "${bins-}" ];
      then
 
			  ${PATH_TO_RELATE}/bin/RelateCoalescentRate \
				  --mode CoalRateForTree \
					--years_per_gen ${years_per_gen} \
          --coal ${output}.coal \
          --dist ${output}.dist \
					-i ${output} \
					-o ${output}

			else

				${PATH_TO_RELATE}/bin/RelateCoalescentRate \
					--mode CoalRateForTree \
					--years_per_gen ${years_per_gen} \
					--bins ${bins} \
          --coal ${output}.coal \
          --dist ${output}.dist \
					-i ${output} \
					-o ${output}

			fi

    done

    if [ -z "${bins-}" ];
    then
     ${PATH_TO_RELATE}/bin/RelateCoalescentRate \
        --mode EstimatePopulationSize \
        --years_per_gen ${years_per_gen} \
        --dist ${output}.dist \
        -i ${output} \
        -o ${output}.pairwise
    else
      ${PATH_TO_RELATE}/bin/RelateCoalescentRate \
        --mode EstimatePopulationSize \
        --years_per_gen ${years_per_gen} \
        --bins ${bins} \
        --dist ${output}.dist \
        -i ${output} \
        -o ${output}.pairwise
    fi

    ${PATH_TO_RELATE}/bin/RelateMutationRate \
      --mode Avg \
      --dist ${output}.dist \
      --years_per_gen ${years_per_gen} \
      -i ${output} \
      -o ${output}

		if [ -z ${noanc-} ]
		then

			#estimate mutation rate and coalescent rates for a final time
			${PATH_TO_RELATE}/bin/RelateExtract \
				--mode DivideAncMut \
				--threads $maxjobs \
				--anc ${filename}.anc \
				--mut ${filename}.mut \
				-o ${output}_tmp

			first_chunk=0
			last_chunk=$(ls ${output}_tmp_chr*.mut.gz | wc -l)
			last_chunk=$((${last_chunk} - 1)) 

			ReEstimateBranchLengths (){

				chunk=$1 
				if [ -z "${seed-}" ];
				then

					${PATH_TO_RELATE}/bin/RelateCoalescentRate \
						--mode ReEstimateBranchLengths \
						--coal ${output}.coal \
						--dist ${output}.dist \
						-m ${mu} \
						-i ${output}_tmp_chr${chunk} \
						-o ${output}_tmp_chr${chunk} 2> ${output}_tmp_chr${chunk}.log 

				else

					${PATH_TO_RELATE}/bin/RelateCoalescentRate \
						--mode ReEstimateBranchLengths \
						--coal ${output}.coal \
						--dist ${output}.dist \
						--seed $(( $seed + ${num_iter} + 1 )) \
						-m ${mu} \
						-i ${output}_tmp_chr${chunk} \
						-o ${output}_tmp_chr${chunk} 2> ${output}_tmp_chr${chunk}.log 

				fi

				rm ${output}_tmp_chr${chunk}.anc.gz
				rm ${output}_tmp_chr${chunk}.mut.gz

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

			parallelize_estimating_branchlengths `seq ${first_chunk} 1 ${last_chunk}`

			for chr in `seq ${first_chunk} 1 ${last_chunk}`
			do
				cat ${output}_tmp_chr${chr}.log  
				rm ${output}_tmp_chr${chr}.log  
			done

			${PATH_TO_RELATE}/bin/RelateExtract \
				--mode CombineAncMut \
				-o ${output}_tmp
				mv ${output}_tmp.anc.gz ${output}.anc.gz
				mv ${output}_tmp.mut.gz ${output}.mut.gz
		else
      rm ${output}.anc
			rm ${output}.mut
			rm ${output}.dist
		fi
  
	else

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

		if [ ! -z "${pop_of_interest-}" ];
		then

		  labels=$(echo ${pop_of_interest} | tr -d ",")

			ExtractSubTrees (){
					chr=$1
					#delete all trees that have fewer than $threshold mutations
					${PATH_TO_RELATE}/bin/RelateExtract \
						--mode SubTreesForSubpopulation \
						--poplabels ${filename_poplabels} \
						--pop_of_interest ${pop_of_interest} \
						--anc ${filename}_chr${chr}.anc \
						--mut ${filename}_chr${chr}.mut \
						-o ${output}_${labels}_chr${chr} 2> ${output}_chr${chr}.log  
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
	  	parallelize_extract_subtrees ${chromosomes}

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


			RemoveTreesWithFewMutations (){
				chr=$1
				#delete all trees that have fewer than $threshold mutations
				${PATH_TO_RELATE}/bin/RelateExtract \
					--mode RemoveTreesWithFewMutations \
					--threshold $threshold \
					--anc ${filename}_chr${chr}.anc \
					--mut ${filename}_chr${chr}.mut \
					-o ${output}_chr${chr} 2> ${output}_chr${chr}.log

				gzip ${output}_chr${chr}.anc
				gzip ${output}_chr${chr}.mut
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
			parallelize_remove_trees ${chromosomes}

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

      if [ -z "${bins-}" ];
      then

        ${PATH_TO_RELATE}/bin/RelateCoalescentRate \
          --mode CoalRateForTree \
          --years_per_gen ${years_per_gen} \
          --dist ${output} \
          --chr ${filename_chr} \
          -i ${output} \
          -o ${output}

      else

        ${PATH_TO_RELATE}/bin/RelateCoalescentRate \
          --mode CoalRateForTree \
          --years_per_gen ${years_per_gen} \
          --bins ${bins} \
          --dist ${output} \
          --chr ${filename_chr} \
          -i ${output} \
          -o ${output}

      fi
	
			#repeat iterations of estimating mutation rate, coalescence rates and re-estimating branch lengths
			for i in `seq 1 1 ${num_iter}`
			do

				for chr in ${chromosomes}
				do

					${PATH_TO_RELATE}/bin/RelateExtract \
						--mode DivideAncMut \
						--threads $maxjobs \
						--anc ${output}_chr${chr}.anc \
						--mut ${output}_chr${chr}.mut \
						-o ${output}_chr${chr}_tmp 

					first_chunk=0
					last_chunk=$(( $(ls ${output}_chr${chr}_tmp_chr*.mut.gz | wc -l) - 1))

					ReEstimateBranchLengths (){

						chunk=$1 
						if [ -z "${seed-}" ];
						then

							${PATH_TO_RELATE}/bin/RelateCoalescentRate \
								--mode SampleBranchLengths \
								--coal ${output}.coal \
								--dist ${output}_chr${chr}.dist \
								--num_samples 1 \
								-m ${mu} \
								-i ${output}_chr${chr}_tmp_chr${chunk} \
								-o ${output}_chr${chr}_tmp_chr${chunk} 2> ${output}_chr${chr}_chr${chunk}.log 

						else

							${PATH_TO_RELATE}/bin/RelateCoalescentRate \
								--mode SampleBranchLengths \
								--coal ${output}.coal \
								--dist ${output}_chr${chr}.dist \
								--num_samples 1 \
								--seed $(( $seed + ${i} )) \
								-m ${mu} \
								-i ${output}_chr${chr}_tmp_chr${chunk} \
								-o ${output}_chr${chr}_tmp_chr${chunk} 2> ${output}_chr${chr}_chr${chunk}.log 
						fi

						rm ${output}_chr${chr}_tmp_chr${chunk}.anc.gz
						rm ${output}_chr${chr}_tmp_chr${chunk}.mut.gz

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
          parallelize_estimating_branchlengths `seq ${first_chunk} 1 ${last_chunk}`

					for chunk in `seq ${first_chunk} 1 ${last_chunk}`
					do
						cat ${output}_chr${chr}_chr${chunk}.log  
						rm ${output}_chr${chr}_chr${chunk}.log  
					done

					${PATH_TO_RELATE}/bin/RelateExtract \
						--mode CombineAncMut \
						-o ${output}_chr${chr}_tmp

					mv ${output}_chr${chr}_tmp.anc.gz ${output}_chr${chr}.anc.gz
					mv ${output}_chr${chr}_tmp.mut.gz ${output}_chr${chr}.mut.gz

	      done


				if [ -z "${bins-}" ];
				then

          ${PATH_TO_RELATE}/bin/RelateCoalescentRate \
            --mode CoalRateForTree \
            --years_per_gen ${years_per_gen} \
            --chr ${filename_chr} \
            --coal ${output}.coal \
            --dist ${output} \
            -i ${output} \
					  -o ${output}	

        else

					${PATH_TO_RELATE}/bin/RelateCoalescentRate \
						--mode CoalRateForTree \
						--years_per_gen ${years_per_gen} \
						--chr ${filename_chr} \
			  		--bins ${bins} \
            --coal ${output}.coal \
            --dist ${output} \
						-i ${output} \
						-o ${output}

				fi

			done

			EstimatePopulationSize (){
				chr=$1
				#delete all trees that have fewer than $threshold mutations
				if [ -z ${bins-} ]
				then
					${PATH_TO_RELATE}/bin/RelateCoalescentRate \
						--mode CoalescentRateForSection \
						--years_per_gen ${years_per_gen} \
            --dist ${output}_chr${chr}.dist \
						-i ${output}_chr${chr} \
						-o ${output}.pairwise_chr${chr} 2> ${output}_chr${chr}.log  
				else
					${PATH_TO_RELATE}/bin/RelateCoalescentRate \
						--mode CoalescentRateForSection \
						--years_per_gen ${years_per_gen} \
						--bins ${bins} \
            --dist ${output}_chr${chr}.dist \
						-i ${output}_chr${chr} \
						-o ${output}.pairwise_chr${chr} 2> ${output}_chr${chr}.log  
				fi
			}

			parallelize_population_size () {
				while [ $# -gt 0 ] ; do
					jobcnt=(`jobs -p`)
					if [ ${#jobcnt[@]} -lt $maxjobs ] ; then
						EstimatePopulationSize $1 &
						shift
					fi
				done
				wait
			}
			parallelize_population_size ${chromosomes}

			for chr in ${chromosomes}
			do
				cat ${output}_chr${chr}.log  
				rm ${output}_chr${chr}.log  
			done

			${PATH_TO_RELATE}/bin/RelateCoalescentRate \
				--mode SummarizeCoalescentRateForGenome \
				--chr ${filename_chr} \
				-o ${output}.pairwise

			${PATH_TO_RELATE}/bin/RelateCoalescentRate \
				--mode FinalizePopulationSize \
				-o ${output}.pairwise

      ${PATH_TO_RELATE}/bin/RelateMutationRate \
        --mode Avg \
        --dist ${output} \
        --years_per_gen ${years_per_gen} \
        --chr ${filename_chr} \
        -i ${output} \
        -o ${output}

		if [ -z ${noanc-} ]
		then

			#estimate mutation rate and coalescent rates for a final time
			for chr in ${chromosomes}
			do

				${PATH_TO_RELATE}/bin/RelateExtract \
					--mode DivideAncMut \
					--threads $maxjobs \
					--anc ${filename}_chr${chr}.anc \
					--mut ${filename}_chr${chr}.mut \
					-o ${output}_chr${chr}_tmp 

				first_chunk=0
				last_chunk=$(( $(ls ${output}_chr${chr}_tmp_chr*.mut.gz | wc -l) - 1))

				ReEstimateBranchLengths (){

					chunk=$1 
					if [ -z "${seed-}" ];
					then

						${PATH_TO_RELATE}/bin/RelateCoalescentRate \
							--mode ReEstimateBranchLengths \
							--coal ${output}.coal \
							--dist ${output}_chr${chr}.dist \
							-m ${mu} \
							-i ${output}_chr${chr}_tmp_chr${chunk} \
							-o ${output}_chr${chr}_tmp_chr${chunk} 2> ${output}_chr${chr}_chr${chunk}.log 

					else

						${PATH_TO_RELATE}/bin/RelateCoalescentRate \
							--mode ReEstimateBranchLengths \
							--coal ${output}.coal \
							--dist ${output}_chr${chr}.dist \
							--seed $(( $seed + ${num_iter} + 1 )) \
							-m ${mu} \
							-i ${output}_chr${chr}_tmp_chr${chunk} \
							-o ${output}_chr${chr}_tmp_chr${chunk} 2> ${output}_chr${chr}_chr${chunk}.log 

					fi

					rm ${output}_chr${chr}_tmp_chr${chunk}.anc.gz
					rm ${output}_chr${chr}_tmp_chr${chunk}.mut.gz

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
				parallelize_estimating_branchlengths `seq ${first_chunk} 1 ${last_chunk}`

				for chunk in `seq ${first_chunk} 1 ${last_chunk}`
				do
					cat ${output}_chr${chr}_chr${chunk}.log  
					rm ${output}_chr${chr}_chr${chunk}.log  
				done

				${PATH_TO_RELATE}/bin/RelateExtract \
					--mode CombineAncMut \
					-o ${output}_chr${chr}_tmp

				mv ${output}_chr${chr}_tmp.anc.gz ${output}_chr${chr}.anc.gz
				mv ${output}_chr${chr}_tmp.mut.gz ${output}_chr${chr}.mut.gz

			done

		else

			for chr in ${chromosomes}
			do
        rm ${output}_chr${chr}.anc.gz
				rm ${output}_chr${chr}.mut.gz
				rm ${output}_chr${chr}.dist
			done

		fi

  fi

fi

if [ ! -z ${filename_poplabels} ];
then
	${PATH_TO_RELATE}/bin/RelateCoalescentRate \
		--mode FinalizePopulationSize \
		-o ${output}.pairwise \
		--poplabels ${filename_poplabels}
fi

if [ ! -z "${pop_of_interest-}" ];
then
  rm ${output}_${labels}*
fi

#plot results
Rscript ${DIR}/plot_population_size.R ${output} ${years_per_gen}


