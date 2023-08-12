#!/bin/bash

#Treat expanding unset parameters as error 
set -u
#Exit with status 1 if a command throws an error
set -e

if [ $# -le 0 ]
then
  echo ""
  echo "Not enough arguments supplied. Execute as"
  echo "PATH_TO_RELATE/scripts/SampleBranchLengths/SampleBranchLengths.sh"
  echo ""
  echo "-i,--input:    Filename of .anc and .mut files without file extension. Needs to contain all SNPs in haps/sample input file." 
  echo "-o, --output:  Filename of output files without file extension."
  echo "-m,--mu:       Mutation rate, e.g., 1.25e-8."
  echo "--coal:        Estimate of coalescent rate, treating all haplotypes as one population."
  echo "--first_bp:        Optional: First bp position."
  echo "--last_bp:         Optional: Last bp position."
  echo "--dist:            Optional: File containing bp and dist between SNPs. Obtained e.g., from RelateExtract --mode AncMutForSubregion. Required if anc/mut file is for subregion."
	echo "--seed:            Optional: Random seed for branch lengths estimation."
	echo "--threads:         Optional. Maximum number of threads."
  echo ""
  exit 1;
fi

PATH_TO_RELATE=$0
PATH_TO_RELATE=$(echo ${PATH_TO_RELATE} | awk -F\scripts/SampleBranchLengths/SampleBranchLengthsParallel.sh '{print $1}')

######################################################################################################

######################## Read arguments from command line ###########################

maxjobs=1
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
    --dist)
      dist="$2"
      shift # past argument
      shift # past value
      ;;
		--threads)
			maxjobs="$2"
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

echo "********************************"
echo "Parameters passed to script:"
echo "input         = $filename"
echo "mu            = $mu"
echo "coal          = $coal" 
echo "output        = $output"

if [ ! -z "${dist-}" ];
then
  echo "dist          = ${dist}"  
fi

if [ ! -z "${first_bp-}" -a ! -z "${last_bp-}" ];
then
  echo "first_bp      = $first_bp"
  echo "last_bp       = $last_bp"
fi

if [ ! -z "${threads-}" ];
then
	echo "threads       = $threads"
fi

if [ ! -z "${seed-}" ];
then
  echo "seed          = $seed"
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

if [ "$maxjobs" -eq 1 ]
then

	if [ ! -z "${first_bp-}" -a ! -z "${last_bp-}" ];
	then

		if [ "${output}" == "${filename}" ]
		then
			echo "Output filename should be different to input filename."
			exit 1
		fi 
		if [ ! -z "${dist-}" ]
		then
			if [ "${output}.dist" == "${dist}" ]
			then
				echo "Output filename should be different to dist filename."
				exit 1
			fi 
		fi

		#create output.anc
		${PATH_TO_RELATE}/bin/RelateExtract \
			--mode AncMutForSubregion \
			--first_bp ${first_bp} \
			--last_bp ${last_bp} \
			--anc ${filename}.anc \
			--mut ${filename}.mut \
			-o ${output}

		if [ ! -z "${dist-}" ]
		then
			is_gzipped=$(file ${dist} | grep -c gzip || true)
			if [ ${is_gzipped} -eq "0" ]
			then
				cp ${dist} ${output}.dist
			else
				gunzip -c ${dist} > ${output}.dist
			fi
		fi

		# sample branch lengths for trees in this region
			
			if [ ! -z "${seed-}" ];
			then

				${PATH_TO_RELATE}/bin/RelateCoalescentRate \
					--mode SampleBranchLengths \
					-m ${mu} \
					--coal ${coal} \
					-i ${output} \
					-o ${output} \
          --num_samples 1 \
					--dist ${output}.dist \
					--seed ${seed}

			else

				${PATH_TO_RELATE}/bin/RelateCoalescentRate \
					--mode SampleBranchLengths \
					-m ${mu} \
					--coal ${coal} \
					-i ${output} \
					-o ${output} \
          --num_samples 1 \
					--dist ${output}.dist 

			fi

	elif [ ! -z "${dist-}" ]
	then

		# sample branch lengths for trees in this region
			if [ ! -z "${seed-}" ];
			then

				${PATH_TO_RELATE}/bin/RelateCoalescentRate \
					--mode SampleBranchLengths \
					-m ${mu} \
					--coal ${coal} \
					--dist ${dist} \
					-i ${filename} \
					-o ${output} \
          --num_samples 1 \
					--seed ${seed}

			else

				${PATH_TO_RELATE}/bin/RelateCoalescentRate \
					--mode SampleBranchLengths \
					-m ${mu} \
					--coal ${coal} \
					--dist ${dist} \
          --num_samples 1 \
					-i ${filename} \
					-o ${output} 

			fi

	else

		# sample branch lengths for trees in this region
			if [ ! -z "${seed-}" ];
			then

				${PATH_TO_RELATE}/bin/RelateCoalescentRate \
					--mode SampleBranchLengths \
					-m ${mu} \
					--coal ${coal} \
          --num_samples 1 \
					-i ${filename} \
					-o ${output} \
					--seed ${seed}

			else

				${PATH_TO_RELATE}/bin/RelateCoalescentRate \
					--mode SampleBranchLengths \
					-m ${mu} \
					--coal ${coal} \
          --num_samples 1 \
					-i ${filename} \
					-o ${output} 

			fi

	fi

else

	if [ ! -z "${first_bp-}" -a ! -z "${last_bp-}" ];
	then

		if [ "${output}" == "${filename}" ]
		then
			echo "Output filename should be different to input filename."
			exit 1
		fi 
		if [ ! -z "${dist-}" ]
		then
			if [ "${output}.dist" == "${dist}" ]
			then
				echo "Output filename should be different to dist filename."
				exit 1
			fi 
		fi

		#create output.anc
		${PATH_TO_RELATE}/bin/RelateExtract \
			--mode AncMutForSubregion \
			--first_bp ${first_bp} \
			--last_bp ${last_bp} \
			--anc ${filename}.anc \
			--mut ${filename}.mut \
			-o ${output}

		if [ ! -z "${dist-}" ]
		then
			is_gzipped=$(file ${dist} | grep -c gzip || true)
			if [ ${is_gzipped} -eq "0" ]
			then
				cp ${dist} ${output}.dist
			else
				gunzip -c ${dist} > ${output}.dist
			fi
		fi

		# sample branch lengths for trees in this region
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
			if [ ! -z "${seed-}" ];
			then

				${PATH_TO_RELATE}/bin/RelateCoalescentRate \
					--mode SampleBranchLengths \
					-m ${mu} \
					--coal ${coal} \
          --num_samples 1 \
					-i ${output}_tmp_chr${chunk} \
					-o ${output}_tmp_chr${chunk} \
					--dist ${output}.dist \
					--seed ${seed} 2> ${output}_tmp_chr${chunk}.log

			else

				${PATH_TO_RELATE}/bin/RelateCoalescentRate \
					--mode SampleBranchLengths \
					-m ${mu} \
					--coal ${coal} \
          --num_samples 1 \
					-i ${output}_tmp_chr${chunk} \
					-o ${output}_tmp_chr${chunk} \
					--dist ${output}.dist 2> ${output}_tmp_chr${chunk}.log

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

	elif [ ! -z "${dist-}" ]
	then

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
			if [ ! -z "${seed-}" ];
			then

				${PATH_TO_RELATE}/bin/RelateCoalescentRate \
					--mode SampleBranchLengths \
					-m ${mu} \
					--coal ${coal} \
          --num_samples 1 \
					-i ${output}_tmp_chr${chunk} \
					-o ${output}_tmp_chr${chunk} \
					--dist ${dist} \
					--seed ${seed} 2> ${output}_tmp_chr${chunk}.log

			else

				${PATH_TO_RELATE}/bin/RelateCoalescentRate \
					--mode SampleBranchLengths \
					-m ${mu} \
					--coal ${coal} \
          --num_samples 1 \
					-i ${output}_tmp_chr${chunk} \
					-o ${output}_tmp_chr${chunk} \
					--dist ${dist} 2> ${output}_tmp_chr${chunk}.log

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

    ${PATH_TO_RELATE}/bin/RelateExtract \
			--mode ExtractDistFromMut \
			--mut ${filename}.mut \
			-o ${output}

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
			if [ ! -z "${seed-}" ];
			then

				${PATH_TO_RELATE}/bin/RelateCoalescentRate \
					--mode SampleBranchLengths \
					-m ${mu} \
					--coal ${coal} \
          --num_samples 1 \
					-i ${output}_tmp_chr${chunk} \
					-o ${output}_tmp_chr${chunk} \
					--dist ${output}.dist \ 
					--seed ${seed} 2> ${output}_tmp_chr${chunk}.log

			else

				${PATH_TO_RELATE}/bin/RelateCoalescentRate \
					--mode SampleBranchLengths \
					-m ${mu} \
					--coal ${coal} \
          --num_samples 1 \
					--dist ${output}.dist \
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

	fi

fi
