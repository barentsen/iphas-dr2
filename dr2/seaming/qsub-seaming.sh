#!/bin/bash

# Crowded regions need more memory to run
for LON in `seq 25 5 115`
do
    qsub -N seaming$LON -v STRIP=$LON -l pmem=4gb pbs-seaming.sh
done

# Less crowded regions
for LON in `seq 120 5 215`
do
    qsub -N seaming$LON -v STRIP=$LON -l pmem=2gb pbs-seaming.sh
done
