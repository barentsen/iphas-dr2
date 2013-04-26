#!/bin/bash
for LON in `seq 30 10 210`
do
    qsub -N seaming$LON -v STRIP=$LON pbs-seaming.sh
done