#!/bin/bash
for BAND in r i ha; do
    qsub -N iphas-offsets-$BAND -v BAND=$BAND pbs-offsets.sh
done
