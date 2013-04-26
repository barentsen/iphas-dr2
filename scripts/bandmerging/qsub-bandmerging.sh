#!/bin/bash
qsub -N bandmerging1 -v LON1=20,LON2=50 pbs-bandmerging.sh
qsub -N bandmerging2 -v LON1=50,LON2=80 pbs-bandmerging.sh
qsub -N bandmerging3 -v LON1=80,LON2=110 pbs-bandmerging.sh
qsub -N bandmerging4 -v LON1=110,LON2=140 pbs-bandmerging.sh
qsub -N bandmerging5 -v LON1=140,LON2=180 pbs-bandmerging.sh
qsub -N bandmerging6 -v LON1=180,LON2=220 pbs-bandmerging.sh
