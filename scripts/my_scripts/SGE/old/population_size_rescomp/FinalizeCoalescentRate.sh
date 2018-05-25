#!/bin/bash

#$ -V
#$ -j y
#$ -P myers.prjc -q short.qc

echo "***********************************************"
echo "SGE job ID: "$JOB_ID
echo "SGE task ID: "$SGE_TASK_ID
echo "Run on host: "`hostname`
echo "Operating system: "`uname -s`
echo "Username: "`whoami`
echo "Started at: "`date`
echo "***********************************************"

../bin/RelateCoalescentRate --mode SummarizeCoalescentRateForGenome --first_chr $first_chr --last_chr $last_chr -o relate_${pop} 
../bin/RelateCoalescentRate --mode FinalizePopulationSize -i relate_${pop} -o relate_${pop} 
../bin/RelateCoalescentRate --mode FinalizePopulationSize --samples ../1000GP_Phase3_${pop}.sample -i relate_${pop} -o relate_${pop}_bypop 
../bin/RelateCoalescentRate --mode FinalizePopulationSize --samples hap -i relate_${pop} -o relate_${pop}_byhap

../bin/RelateMutationRate --mode SummarizeForGenome --first_chr $first_chr --last_chr $last_chr -o relate_${pop} 
../bin/RelateMutationRate --mode Finalize -i relate_${pop} -o relate_${pop} 
../bin/RelateMutationRate --mode FinalizeAvg -i relate_${pop} -o relate_${pop}_avg

../bin/RelateMutationRate --mode Avg --pos pos_chr1_${pop}.txt --mut relate_chr1_${pop}_pre.mut -i relate_chr1_${pop} -o relate_chr1_${pop}_avg

module load R/3.4.0
Rscript ../update_coal.R ${pop}
Rscript ../plot_population_size.R ${pop}
Rscript ../MutationRateOverTime.R ${pop}


echo "***********************************************"
echo "Finished at: "`date`
echo "***********************************************"
exit 0

