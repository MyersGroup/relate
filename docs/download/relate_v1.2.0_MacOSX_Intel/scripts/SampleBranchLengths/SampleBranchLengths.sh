#!/bin/bash

#Treat expanding unset parameters as error 
set -u
#Exit with status 1 if a command throws an error
set -e

if [ $# -le 0 ]
then
  echo ""
  echo "Not enough arguments supplied. Execute as"
  echo "PATH_TO_RELATE/scripts/SampleBrancheLengths/SampleBranchLengths.sh"
  echo ""
  echo "-i,--input:    Filename of .anc and .mut files without file extension. Needs to contain all SNPs in haps/sample input file." 
  echo "-o, --output:  Filename of output files without file extension."
  echo "-m,--mu:       Mutation rate, e.g., 1.25e-8."
  echo "--coal:        Estimate of coalescent rate, treating all haplotypes as one population."
  echo "--num_samples: Number of times branch lengths are sampled."
  echo "--first_bp:        Optional: First bp position."
  echo "--last_bp:         Optional: Last bp position."
  echo "--dist:            Optional: File containing bp and dist between SNPs. Obtained e.g., from RelateExtract --mode AncMutForSubregion. Required if anc/mut file is for subregion."
  echo "--num_proposals:   Optional: Number of MCMC proposals between samples. Default is max(10*N, 1000), where N is the number of haplotypes. If num_proposals=0, all samples will be identical to the input."
	echo "--format:          Optional: Output file format when sampling branches. a: anc/mut, n: newick, b: binary. Default a."
	echo "--seed:            Optional: Random seed for branch lengths estimation."
  echo ""
  exit 1;
fi

PATH_TO_RELATE=$0
PATH_TO_RELATE=$(echo ${PATH_TO_RELATE} | awk -F\scripts/SampleBranchLengths/SampleBranchLengths.sh '{print $1}')

######################################################################################################

######################## Read arguments from command line ###########################

format=a

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
    --num_samples)
      num_samples="$2"
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
    --num_proposals)
      num_proposals="$2"
      shift # past argument
      shift # past value
      ;;
		--format)
			format="$2"
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
echo "num_samples   = $num_samples" 
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

if [ ! -z "${num_proposals-}" ];
then
  echo "num_proposals = ${num_proposals}"  
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
  # this generates .newick, .sites
  if [ ! -z "${num_proposals-}" ];
  then
    if [ ! -z "${seed-}" ];
    then

      ${PATH_TO_RELATE}/bin/RelateCoalescentRate \
        --mode SampleBranchLengths \
        -m ${mu} \
        --coal ${coal} \
        --num_samples ${num_samples} \
        -i ${outp2ut} \
        -o ${output} \
        --num_proposals ${num_proposals} \
        --dist ${output}.dist \
				--format ${format} \
        --seed ${seed}

    else

      ${PATH_TO_RELATE}/bin/RelateCoalescentRate \
        --mode SampleBranchLengths \
        -m ${mu} \
        --coal ${coal} \
        --num_samples ${num_samples} \
        --num_proposals ${num_proposals} \
        -i ${output} \
        -o ${output} \
				--format ${format} \
        --dist ${output}.dist 

    fi
  else
    if [ ! -z "${seed-}" ];
    then

      ${PATH_TO_RELATE}/bin/RelateCoalescentRate \
        --mode SampleBranchLengths \
        -m ${mu} \
        --coal ${coal} \
        --num_samples ${num_samples} \
        -i ${output} \
        -o ${output} \
        --dist ${output}.dist \
				--format ${format} \
        --seed ${seed}

    else

      ${PATH_TO_RELATE}/bin/RelateCoalescentRate \
        --mode SampleBranchLengths \
        -m ${mu} \
        --coal ${coal} \
        --num_samples ${num_samples} \
        -i ${output} \
        -o ${output} \
				--format ${format} \
        --dist ${output}.dist 

    fi
  fi

elif [ ! -z "${dist-}" ]
then

  # sample branch lengths for trees in this region
  # this generates .newick, .sites
  if [ ! -z "${num_proposals-}" ];
  then
    if [ ! -z "${seed-}" ];
    then

      ${PATH_TO_RELATE}/bin/RelateCoalescentRate \
        --mode SampleBranchLengths \
        -m ${mu} \
        --coal ${coal} \
        --num_samples ${num_samples} \
        --num_proposals ${num_proposals} \
        --dist ${dist} \
        -i ${filename} \
        -o ${output} \
				--format ${format} \
        --seed ${seed}

    else

      ${PATH_TO_RELATE}/bin/RelateCoalescentRate \
        --mode SampleBranchLengths \
        -m ${mu} \
        --coal ${coal} \
        --num_samples ${num_samples} \
        --num_proposals ${num_proposals} \
        --dist ${dist} \
				--format ${format} \
        -i ${filename} \
        -o ${output} 

    fi
  else
    if [ ! -z "${seed-}" ];
    then

      ${PATH_TO_RELATE}/bin/RelateCoalescentRate \
        --mode SampleBranchLengths \
        -m ${mu} \
        --coal ${coal} \
        --num_samples ${num_samples} \
        --dist ${dist} \
				--format ${format} \
        -i ${filename} \
        -o ${output} \
        --seed ${seed}

    else

      ${PATH_TO_RELATE}/bin/RelateCoalescentRate \
        --mode SampleBranchLengths \
        -m ${mu} \
        --coal ${coal} \
        --num_samples ${num_samples} \
        --dist ${dist} \
				--format ${format} \
        -i ${filename} \
        -o ${output} 

    fi
  fi

else

  # sample branch lengths for trees in this region
  # this generates .newick, .sites
  if [ ! -z "${num_proposals-}" ];
  then
    if [ ! -z "${seed-}" ];
    then

      ${PATH_TO_RELATE}/bin/RelateCoalescentRate \
        --mode SampleBranchLengths \
        -m ${mu} \
        --coal ${coal} \
        --num_samples ${num_samples} \
        --num_proposals ${num_proposals} \
				--format ${format} \
        -i ${filename} \
        -o ${output} \
        --seed ${seed}

    else

      ${PATH_TO_RELATE}/bin/RelateCoalescentRate \
        --mode SampleBranchLengths \
        -m ${mu} \
        --coal ${coal} \
        --num_samples ${num_samples} \
        --num_proposals ${num_proposals} \
				--format ${format} \
        -i ${filename} \
        -o ${output} 

    fi
  else
    if [ ! -z "${seed-}" ];
    then

      ${PATH_TO_RELATE}/bin/RelateCoalescentRate \
        --mode SampleBranchLengths \
        -m ${mu} \
        --coal ${coal} \
        --num_samples ${num_samples} \
				--format ${format} \
        -i ${filename} \
        -o ${output} \
        --seed ${seed}

    else

      ${PATH_TO_RELATE}/bin/RelateCoalescentRate \
        --mode SampleBranchLengths \
        -m ${mu} \
        --coal ${coal} \
        --num_samples ${num_samples} \
				--format ${format} \
        -i ${filename} \
        -o ${output} 

    fi
  fi

fi

