#!/bin/bash
# Executes the band-merging script
# distributed by splitting up the work in 10-deg longitude strips
for LON1 in `seq 20 10 210`
do
    LON2=$(expr $LON1 + 10)
    qsub -N bandmerging$LON1 -v LON1=$LON1,LON2=$LON2 pbs-bandmerging.sh
done

#qsub -N bandmerging1 -v LON1=20,LON2=40 pbs-bandmerging.sh

