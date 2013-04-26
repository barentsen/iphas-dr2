#!/bin/bash
#qsub -N seaming -v STRIP=$LON pbs-seaming.sh

for LON in `seq 30 10 210`
do
    echo $LON
done