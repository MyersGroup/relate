#!/bin/bash

#Treat expanding unset parameters as error 
set -u
#Exit with status 1 if a command throws an error
set -e

if [ $# -le 0 ]
then
  echo ""
  echo "Not enough arguments supplied. Execute as"
  echo "PATH_TO_RELATE/scripts/PrepareInputFiles/PrepareInputFiles.sh"
  echo ""
  echo "--haps:       Filename of haps file."
  echo "--sample:     Filename of sample file."
  echo "--ancestor:   Filename of fasta file containing ancestral genome."
  echo "-o, --output: Filename of output files without file extension."
  echo "--mask:       Optional: Filename of fasta file of same length as --ancestor containing genome mask."
  echo "--remove_ids: Optional: Filename of file containing sample ids to be removed from data. (One id per line)"
  echo "--poplabels:  Optional: Filename of file containing population labels. Samples must be listed in the same order as the .sample file."
  echo ""
  exit 1;
fi

PATH_TO_RELATE=$0
PATH_TO_RELATE=$(echo ${PATH_TO_RELATE} | awk -F\scripts/PrepareInputFiles/PrepareInputFiles.sh '{print $1}')

######################################################################################################

######################## Read arguments from command line ###########################

maxjobs=1

POSITIONAL=()
while [[ $# -gt 0 ]]
do
  key="$1"

  case $key in
    --haps)
      filename_haps="$2"
      shift # past argument
      shift # past value
      ;;
    --sample)
      filename_sample="$2"
      shift # past argument
      shift # past value
      ;;
    --ancestor)
      filename_ancestor="$2"
      shift # past argument
      shift # past value
      ;;
    --remove_ids)
      filename_remove="$2"
      shift # past argument
      shift # past value
      ;;
    --mask)
      filename_mask="$2"
      shift # past argument
      shift # past value
      ;;
    --poplabels)
      filename_poplabels="$2"
      shift # past argument
      shift # past value
      ;;
    -o|--output)
      filename_output="$2"
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
echo "haps        = ${filename_haps}"
echo "sample      = ${filename_sample}"
echo "ancestor    = ${filename_ancestor}"
echo "output      = ${filename_output}"
if [ ! -z "${filename_mask-}" ];
then
  echo "mask        = $filename_mask"
fi
if [ ! -z "${filename_remove-}" ];
then
  echo "remove_ids  = $filename_remove"
fi
if [ ! -z "${filename_poplabels-}" ];
then
  echo "poplabels   = $filename_poplabels"
fi
echo "********************************"


# remove biallelic SNPs
${PATH_TO_RELATE}/bin/RelateFileFormats \
  --mode RemoveNonBiallelicSNPs \
  --haps ${filename_haps} \
  -o ${filename_output}_biall

# flip haps using ancestor
${PATH_TO_RELATE}/bin/RelateFileFormats \
  --mode FlipHapsUsingAncestor \
  --haps ${filename_output}_biall.haps \
  --sample ${filename_sample} \
  -o ${filename_output}_ancest \
  --ancestor ${filename_ancestor}

rm ${filename_output}_biall.haps

is_gzipped=$(file ${filename_sample} | grep -c "gzip" || true)
if [ ${is_gzipped} -eq 0 ];
then
  cp ${filename_sample} ${filename_output}.sample
else
  gunzip -c ${filename_sample} > ${filename_output}.sample
fi

if [ ! -z "${filename_remove-}" ];
then

  if [ ! -z "${filename_poplabels-}" ];
  then
    # Subset individuals
    ${PATH_TO_RELATE}/bin/RelateFileFormats \
      --mode RemoveSamples \
      --haps ${filename_output}_ancest.haps \
      --sample ${filename_output}.sample \
      --poplabels ${filename_poplabels} \
      -i ${filename_remove} \
      -o ${filename_output}_rem

    mv ${filename_output}_rem.poplabels ${filename_output}.poplabels
    filename_poplabels=${filename_output}.poplabels
  else
    # Subset individuals
    ${PATH_TO_RELATE}/bin/RelateFileFormats \
      --mode RemoveSamples \
      --haps ${filename_output}_ancest.haps \
      --sample ${filename_output}.sample \
      -i ${filename_remove} \
      -o ${filename_output}_rem
  fi


  rm ${filename_output}_ancest.haps
  mv ${filename_output}_rem.sample ${filename_output}.sample

  if [ ! -z "${filename_mask-}" ];
  then

    # filter data using genomic mask
    ${PATH_TO_RELATE}/bin/RelateFileFormats \
      --mode FilterHapsUsingMask \
      --haps ${filename_output}_rem.haps \
      --sample ${filename_output}.sample \
      -o ${filename_output}_filtered \
      --mask ${filename_mask}

    rm ${filename_output}_rem.haps
    mv ${filename_output}_filtered.haps ${filename_output}.haps
    mv ${filename_output}_filtered.dist ${filename_output}.dist

  else
    mv ${filename_output}_rem.haps ${filename_output}.haps
  fi

else

  if [ ! -z "${filename_mask-}" ];
  then

    # filter data using genomic mask
    ${PATH_TO_RELATE}/bin/RelateFileFormats \
      --mode FilterHapsUsingMask \
      --haps ${filename_output}_ancest.haps \
      --sample ${filename_output}.sample \
      -o ${filename_output}_filtered \
      --mask ${filename_mask}

    rm ${filename_output}_ancest.haps
    mv ${filename_output}_filtered.haps ${filename_output}.haps
    mv ${filename_output}_filtered.dist ${filename_output}.dist

  else
    mv ${filename_output}_ancest.haps ${filename_output}.haps
  fi

fi

if [ ! -z "${filename_poplabels-}" ];
then
    ${PATH_TO_RELATE}/bin/RelateFileFormats \
      --mode GenerateSNPAnnotations \
      --haps ${filename_output}.haps \
      --sample ${filename_output}.sample \
      --poplabels ${filename_poplabels} \
      --ancestor ${filename_ancestor} \
      -o ${filename_output} 
fi

gzip ${filename_output}.haps
gzip ${filename_output}.sample
gzip ${filename_output}.dist

