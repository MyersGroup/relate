#!/bin/bash

#$ -cwd -V
#$ -j y
#$ -P myers.prjc -q short.qc
#$ -e error.out
#$ -o std.out

echo "***********************************************"
echo "SGE job ID: "$JOB_ID
echo "SGE task ID: "$SGE_TASK_ID
echo "Run on host: "`hostname`
echo "Operating system: "`uname -s`
echo "Username: "`whoami`
echo "Started at: "`date`
echo "***********************************************"

gunzip -c ${path_to_input}${input}_chr$c.arg.gz > ${input}_chr$c.arg
gunzip -c ${path_to_input}${input}_chr$c.mut.gz > ${input}_chr$c.mut 

echo "***********************************************"
echo "Finished at: "`date`
echo "***********************************************"
exit 0

