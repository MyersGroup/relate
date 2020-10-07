#!/bin/bash

##1. Read arguments from command line
##2. Define a function that applies stages 2-6 of the algorithm to each chunk
##3. Apply stage 1
##4. For loop over chunks
##5. Apply stage 7

#Treat expanding unset parameters as error 
set -u
#Exit with status 1 if a command throws an error
set -e


if [ $# -le 0 ]
then
  echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
  echo "Not enough arguments supplied. Execute as"
  echo "PATH_TO_RELATE/scripts/RelateParallel/RelateParallel.sh"
  echo ""
  echo "--haps:   Filename of haps file (Output file format of Shapeit)."
  echo "--sample: Filename of sample file (Output file format of Shapeit)."
  echo "--map:    Genetic map."
  echo "-m,--mu:  Mutation rate."
  echo "-N,--Ne:  Effective population size."
  echo "-o,--output: Filename of output without file extension."
  echo "--dist:   Optional but recommended. Distance in BP between SNPs. Can be generated using RelateFileFormats. If unspecified, distances in haps are used."
  echo "--annot:  Optional. Filename of file containing additional annotation of snps. Can be generated using RelateFileFormats."
  echo "--memory: Optional. Approximate memory allowance in GB for storing distance matrices. Default is 5GB."
  echo "sample_ages: Optional. Filename of file containing sample ages (one per line)." 
  echo "--coal:   Optional. Filename of file containing coalescent rates. If specified, it will overwrite --effectiveN." 
  echo "--painting: Optional. Copying and transition parameters in
                    chromosome painting algorithm. Format: theta,rho. Default: 0.025,1." 
  echo "--seed:   Optional. Seed for MCMC in branch lengths estimation."
  echo "--threads:Optional. Maximum number of threads."
  exit 1;
fi

PATH_TO_RELATE=$0
PATH_TO_RELATE=$(echo ${PATH_TO_RELATE} | awk -F\scripts/RelateParallel/RelateParallel.sh '{print $1}')

######################################################################################################

######################## Read arguments from command line ###########################

maxjobs=1
painting="0.025,1"

POSITIONAL=()
while [[ $# -gt 0 ]]
do
  key="$1"

  case $key in
    --haps)
      haps="$2"
      shift # past argument
      shift # past value
      ;;
    --sample)
      sample="$2"
      shift # past argument
      shift # past value
      ;;
    --map)
      map="$2"
      shift # past argument
      shift # past value
      ;;
    -m|--mu)
      mu="$2"
      shift # past argument
      shift # past value
      ;;
    -N|--Ne)
      Ne="$2"
      shift # past argument
      shift # past value
      ;;
    -o|--output)
      output="$2"
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
    --dist)
      dist="$2"
      shift # past argument
      shift # past value
      ;;
    --annot)
      annot="$2"
      shift # past argument
      shift # past value
      ;;
    --memory)
      memory="$2"
      shift # past argument
      shift # past value
      ;;
    --sample_ages)
      sample_ages="$2"
      shift # past argument
      shift # past value
      ;;
		--painting)
			painting="$2"
			shift # past argument
			shift # past value
			;;
    --coal)
      coal="$2"
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
echo "Parameters passed to Relate:"
echo "haps   = $haps"
echo "sample = $sample"
echo "map    = $map"
if [ ! -z "${coal-}" ];
then
  echo "coal   = $coal"
else 
  echo "Ne     = $Ne"
fi
echo "mu     = $mu"
echo "output = $output"
if [ ! -z "${dist-}" ];
then
  echo "dist   = $dist"
fi
if [ ! -z "${annot-}" ];
then
  echo "annot  = $annot"
fi
if [ ! -z "${memory-}" ];
then
  echo "memory = $memory"
fi
if [ ! -z "${painting-}" ];
then
  echo "painting   = $painting"
fi
if [ ! -z "${sample_ages-}" ];
then
  echo "sample_ages = $sample_ages"
fi
if [ ! -z "${seed-}" ];
then
  echo "seed   = $seed"
fi
echo "Maximum number of threads: $maxjobs" 
echo "********************************"


######################################################################################################

########### Function that takes chunk_index and applies Relate stages 2-6 to it ######################

RelateForChunk (){
  chunk_index=$1

  ## paint all sequences against each other
  ${PATH_TO_RELATE}/bin/Relate \
    --mode Paint \
    --painting ${painting} \
    -o ${output} \
    --chunk_index ${chunk_index} 

  #GET NUMBER OF SECTIONS
  num_sections=$(ls ${output}/chunk_${chunk_index}/paint/*bin | wc -l)

  ## build tree topologies
  parallelize_tree_building () {
    while [ $# -gt 0 ] ; do
      jobcnt=(`jobs -p`)
      if [ ${#jobcnt[@]} -lt $maxjobs ] ; then
				if [ -z "${seed-}" ];
				then
					${PATH_TO_RELATE}/bin/Relate \
						--mode BuildTopology \
						--chunk_index ${chunk_index} \
						--first_section $1 \
						--last_section $1 \
						--painting ${painting} \
						-o ${output} 2> ${output}/chunk_${chunk_index}/sec${1}.log &
				else
					${PATH_TO_RELATE}/bin/Relate \
						--mode BuildTopology \
						--chunk_index ${chunk_index} \
						--first_section $1 \
						--last_section $1 \
						--seed ${seed} \
						--painting ${painting} \
						-o ${output} 2> ${output}/chunk_${chunk_index}/sec${1}.log &
				fi
        shift
      fi
    done
    wait
  }
  parallelize_tree_building `seq 0 1 $(($num_sections - 1))`

  for section in `seq 0 1 $(($num_sections - 1))`
  do
    cat ${output}/chunk_${chunk_index}/sec${section}.log 
    rm ${output}/chunk_${chunk_index}/sec${section}.log 
  done

  ## find equivalent branches in adjacent trees 
  ${PATH_TO_RELATE}/bin/Relate \
    --mode FindEquivalentBranches \
    --chunk_index ${chunk_index} \
    -o ${output} 

  ## infer branch lengths
  parallelize_branch_lengths () {
    while [ $# -gt 0 ] ; do
      jobcnt=(`jobs -p`)
      if [ ${#jobcnt[@]} -lt $maxjobs ] ; then

				if [ -z "${seed-}" ];
				then
          if [ ! -z "${sample_ages-}" ]; then

						if [ ! -z "${coal-}" ]; then
							${PATH_TO_RELATE}/bin/Relate \
								--mode InferBranchLengths \
								-m $mu \
								--coal $coal \
								--chunk_index ${chunk_index} \
								--first_section $1 \
								--last_section $1 \
								--sample_ages ${sample_ages} \
								-o ${output} 2> ${output}/chunk_${chunk_index}/sec${1}.log &
						else
							${PATH_TO_RELATE}/bin/Relate \
								--mode InferBranchLengths \
								-m $mu \
								-N $Ne \
								--chunk_index ${chunk_index} \
								--first_section $1 \
								--last_section $1 \
								--sample_ages ${sample_ages} \
								-o ${output} 2> ${output}/chunk_${chunk_index}/sec${1}.log &
						fi

					else

						if [ ! -z "${coal-}" ]; then
							${PATH_TO_RELATE}/bin/Relate \
								--mode InferBranchLengths \
								-m $mu \
								--coal $coal \
								--chunk_index ${chunk_index} \
								--first_section $1 \
								--last_section $1 \
								-o ${output} 2> ${output}/chunk_${chunk_index}/sec${1}.log &
						else
							${PATH_TO_RELATE}/bin/Relate \
								--mode InferBranchLengths \
								-m $mu \
								-N $Ne \
								--chunk_index ${chunk_index} \
								--first_section $1 \
								--last_section $1 \
								-o ${output} 2> ${output}/chunk_${chunk_index}/sec${1}.log &
						fi

					fi

				else

					if [ ! -z "${sample_ages-}" ]; then

						if [ ! -z "${coal-}" ]; then
							${PATH_TO_RELATE}/bin/Relate \
								--mode InferBranchLengths \
								-m $mu \
								--coal $coal \
								--chunk_index ${chunk_index} \
								--first_section $1 \
								--last_section $1 \
								--sample_ages ${sample_ages} \
								--seed ${seed} \
								-o ${output} 2> ${output}/chunk_${chunk_index}/sec${1}.log &
						else
							${PATH_TO_RELATE}/bin/Relate \
								--mode InferBranchLengths \
								-m $mu \
								-N $Ne \
								--chunk_index ${chunk_index} \
								--first_section $1 \
								--last_section $1 \
								--sample_ages ${sample_ages} \
								--seed ${seed} \
								-o ${output} 2> ${output}/chunk_${chunk_index}/sec${1}.log &
						fi

					else

						if [ ! -z "${coal-}" ]; then
							${PATH_TO_RELATE}/bin/Relate \
								--mode InferBranchLengths \
								-m $mu \
								--coal $coal \
								--chunk_index ${chunk_index} \
								--first_section $1 \
								--last_section $1 \
								--seed ${seed} \
								-o ${output} 2> ${output}/chunk_${chunk_index}/sec${1}.log &
						else
						  ${PATH_TO_RELATE}/bin/Relate \
							  --mode InferBranchLengths \
								-m $mu \
								-N $Ne \
								--chunk_index ${chunk_index} \
								--first_section $1 \
								--last_section $1 \
								--seed ${seed} \
								-o ${output} 2> ${output}/chunk_${chunk_index}/sec${1}.log &
						fi

					fi

        fi

        shift
      fi
    done
    wait
  }
  parallelize_branch_lengths `seq 0 1 $(($num_sections - 1))`

  for section in `seq 0 1 $(($num_sections - 1))`
  do
    cat ${output}/chunk_${chunk_index}/sec${section}.log 
    rm ${output}/chunk_${chunk_index}/sec${section}.log 
  done


  ## combine args into one file
  ${PATH_TO_RELATE}/bin/Relate \
    --mode CombineSections \
    -N $Ne \
    --chunk_index ${chunk_index} \
    -o ${output} 
} 


########################### FUNCTION THAT PARALLELISES CODE #######

parallelize () {
  while [ $# -gt 0 ] ; do
    jobcnt=(`jobs -p`)
    if [ ${#jobcnt[@]} -lt $maxjobs ] ; then
      RelateForChunk $1 &
      shift
    fi
  done
  wait
}



######################################################################################################

############################ START RELATE ########################

########## Divide data into chunks (the chunks are overlapping to avoid edge effects) ################
# If statements for optional arguments passed to Relate

if [ -z "${dist-}" ];
then
  if [ -z "${memory-}" ];
  then

      # no optional arguments set
      ${PATH_TO_RELATE}/bin/Relate \
        --mode "MakeChunks" \
        --haps ${haps} \
        --sample ${sample} \
        -o ${output} \
        --map ${map}  

  else

      # only memory is set
      ${PATH_TO_RELATE}/bin/Relate \
        --mode "MakeChunks" \
        --haps ${haps} \
        --sample ${sample} \
        --map ${map} \
        -o ${output} \
        --memory ${memory}

  fi
else
  if [ -z "${memory-}" ];
  then

      # only dist is set
      ${PATH_TO_RELATE}/bin/Relate \
        --mode "MakeChunks" \
        --haps ${haps} \
        --sample ${sample} \
        --map ${map} \
        -o ${output} \
        --dist ${dist}
 
  else

      # dist and memory are set
      ${PATH_TO_RELATE}/bin/Relate \
        --mode "MakeChunks" \
        --haps ${haps} \
        --sample ${sample} \
        --map ${map} \
        --dist ${dist} \
        -o ${output} \
        --memory ${memory}

  fi
fi

###########GET NUMBER OF CHUNKS
num_chunks=$(ls ${output}/parameters_* | wc -l)

echo "********************************"
echo "Number of chunks: "$num_chunks
echo "********************************"

############################ APPLY CODE TO EACH CHUNK ############

chunks=`seq 0 $(($num_chunks - 1))`
for c in $chunks
do
  RelateForChunk $c 
done
###########
#if you want to parallelize this for loop use this:
#parallelize $chunks

############################# Finalize ###########################

#if statement to include annot if specified

if [ -z "${sample_ages-}" ];
then
  if [ -z "${annot-}" ];
  then

    #annot not set
    ${PATH_TO_RELATE}/bin/Relate \
      --mode "Finalize" \
      -o ${output} 

  else

    #annot is set
    ${PATH_TO_RELATE}/bin/Relate \
      --mode "Finalize" \
      --annot ${annot} \
      -o ${output} 

  fi
else
  if [ -z "${annot-}" ];
  then

    #annot not set
    ${PATH_TO_RELATE}/bin/Relate \
      --mode "Finalize" \
      --sample_ages ${sample_ages} \
      -o ${output} 

  else

    #annot is set
    ${PATH_TO_RELATE}/bin/Relate \
      --mode "Finalize" \
      --annot ${annot} \
      --sample_ages ${sample_ages} \
      -o ${output} 

  fi
fi
