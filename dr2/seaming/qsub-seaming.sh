#!/bin/bash
# Crowded fields need more memory to run
for LON in `seq 30 10 100`
do
    qsub -N seaming$LON -v STRIP=$LON -l pmem=4gb pbs-seaming.sh
done

for LON in `seq 110 10 210`
do
    qsub -N seaming$LON -v STRIP=$LON -l pmem=2gb pbs-seaming.sh
done
