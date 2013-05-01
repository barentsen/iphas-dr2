#!/bin/bash

# Crowded regions need more memory to run
for LON in `seq 25 5 115`
do
    qsub -N seaming$LON -v STRIP=$LON -l pmem=3gb,walltime=08:00:00 pbs-seaming.sh
done

# Less crowded regions
for LON in `seq 120 5 215`
do
    qsub -N seaming$LON -v STRIP=$LON -l pmem=2gb,walltime=06:00:00 pbs-seaming.sh
done
