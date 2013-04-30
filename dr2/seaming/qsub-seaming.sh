#!/bin/bash
# Crowded fields need more memory to run
for LON in `seq 25 5 100`
do
    qsub -N seaming$LON -v STRIP=$LON -l pmem=4gb pbs-seaming.sh
done

for LON in `seq 110 5 215`
do
    qsub -N seaming$LON -v STRIP=$LON -l pmem=2gb pbs-seaming.sh
done
