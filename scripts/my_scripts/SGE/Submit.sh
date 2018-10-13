#!/bin/bash

for c in {3..20}
do
  qsub -cwd \
    -hold_jid PrepareData2 \
    -N PrepareData2 \
    -v c=${c} \
    -o PrepareData.log \
    -e PrepareData.log \
    PrepareData.sh
done

