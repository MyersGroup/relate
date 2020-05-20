#!/bin/bash

#Treat expanding unset parameters as error 
set -u
#Exit with status 1 if a command throws an error
set -e

if [ $# -le 0 ]
then
  echo ""
  echo "Not enough arguments supplied. Execute as"
  echo "PATH_TO_RELATE/scripts/TreeView/TreeView.sh"
  echo ""
  echo "--haps:           Filename of haps file."
  echo "--sample:         Filename of sample file."
  echo "--poplabels:      Optional: Filename of file containing population labels. Samples must be listed in the same order as the .sample file."
  echo "--anc:            Filename of anc file."
  echo "--mut:            Filename of mut file."
  echo "--bp_of_interest: Base pair position at which tree will be plotted."
  echo "--years_per_gen : Years per generation."
  echo "-o, --output:     Filename of plot without file extension."
  echo ""
  exit 1;
fi

PATH_TO_RELATE=$0
PATH_TO_RELATE=$(echo ${PATH_TO_RELATE} | awk -F\scripts/TreeView/TreeView.sh '{print $1}')

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
    --anc)
      filename_anc="$2"
      shift # past argument
      shift # past value
      ;;
    --mut)
      filename_mut="$2"
      shift # past argument
      shift # past value
      ;;
    --bp_of_interest)
      bp_of_interest="$2"
      shift # past argument
      shift # past value
      ;;
    --poplabels)
      filename_poplabels="$2"
      shift # past argument
      shift # past value
      ;;
    --years_per_gen)
      years_per_gen="$2"
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
echo "haps           = ${filename_haps}"
echo "sample         = ${filename_sample}"
echo "poplabels      = $filename_poplabels"
echo "anc            = ${filename_anc}"
echo "mut            = ${filename_mut}"
echo "years_per_gen  = ${years_per_gen}"
echo "bp_of_interest = ${bp_of_interest}"
echo "output         = ${filename_output}.pdf"
echo "********************************"

Rscript ${PATH_TO_RELATE}/scripts/TreeView/treeview.R \
  ${PATH_TO_RELATE} ${filename_haps} ${filename_sample} $filename_poplabels ${filename_anc} ${filename_mut} ${years_per_gen} ${bp_of_interest} ${filename_output}

