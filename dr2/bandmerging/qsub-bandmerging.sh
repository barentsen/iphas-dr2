#!/bin/bash
for LON1 in `seq 20 10 210`
do
    LON2=$(expr $LON1 + 10)
    qsub -N bandmerging$LON1 -v LON1=$LON1,LON2=$LON2 pbs-bandmerging.sh
done

#qsub -N bandmerging1 -v LON1=20,LON2=40 pbs-bandmerging.sh
#qsub -N bandmerging2 -v LON1=40,LON2=60 pbs-bandmerging.sh
#qsub -N bandmerging3 -v LON1=60,LON2=90 pbs-bandmerging.sh
#qsub -N bandmerging4 -v LON1=90,LON2=120 pbs-bandmerging.sh
#qsub -N bandmerging5 -v LON1=120,LON2=170 pbs-bandmerging.sh
#qsub -N bandmerging6 -v LON1=170,LON2=220 pbs-bandmerging.sh
