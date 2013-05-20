#!/bin/sh -f    
#PBS -N calibration                                                       
#PBS -m aeb 
#PBS -M gb
#PSB -l pmem=1gb
#PBS -l nodes=1:ppn=8
#PBS -k oe                                                                      
#PBS -q cmain
#PBS -l walltime=23:59:00                                                       

echo $1

echo ------------------------------------------------------
echo -n 'Job is running on node '; cat $PBS_NODEFILE
echo ------------------------------------------------------
echo PBS: qsub is running on $PBS_O_HOST
echo PBS: originating queue is $PBS_O_QUEUE
echo PBS: executing queue is $PBS_QUEUE
echo PBS: working directory is $PBS_O_WORKDIR
echo PBS: execution mode is $PBS_ENVIRONMENT
echo PBS: job identifier is $PBS_JOBID
echo PBS: job name is $PBS_JOBNAME
echo PBS: node file is $PBS_NODEFILE
echo PBS: array ID is $PBS_ARRAYID
echo PBS: current home directory is $PBS_O_HOME
echo PBS: PATH = $PBS_O_PATH
echo ------------------------------------------------------

export PYTHONPATH=/home/gb/bin/epd-7.3-1-rh5-x86_64/lib/python2.7/site-packages/
export PATH=/home/gb/bin/epd-7.3-1-rh5-x86_64/bin:$PATH
export LD_LIBRARY_PATH=/home/gb/bin/wcslib-4.15/lib:$LD_LIBRARY_PATH
cd /home/gb/dev/iphas-dr2/dr2
nice python calibration.py 
echo ------------------------------------------------------
echo Job ends
