#!/bin/bash

#DATADIR='/media/0133d764-0bfe-4007-a9cc-a7b1f61c4d1d/iphas'

DATADIR='/car-data/gb/iphas'
for subdir in `ls $DATADIR`; do
    echo $subdir
    qsub -v IPHASDIR=$subdir iphasDetection.pbs
done
