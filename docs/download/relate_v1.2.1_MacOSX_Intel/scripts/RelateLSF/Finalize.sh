#!/bin/bash

#$ -V
#$ -j y

echo "***********************************************"
echo "LSF job ID: "$LSB_JOBID
echo "LSF task ID: "$LSB_JOBINDEX
echo "Run on host: "`hostname`
echo "Operating system: "`uname -s`
echo "Username: "`whoami`
echo "Started at: "`date`
echo "***********************************************"

if [ -z ${sample_ages-} ]
then
	${PATH_TO_RELATE}/bin/Relate \
		--mode "Finalize" \
		-o ${output} 2>> log/combine_args.log
else
	${PATH_TO_RELATE}/bin/Relate \
		--mode "Finalize" \
		--sample_ages ${sample_ages} \
		-o ${output} 2>> log/combine_args.log
fi

echo "***********************************************"
echo "Finished at: "`date`
echo "***********************************************"
exit 0

