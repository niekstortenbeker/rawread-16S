#!/bin/bash

#$ -S /bin/bash
#$ -cwd
#$ -l scratch
#$ -V
#$ -pe threaded 8
#$ -R y
#$ -e bash_error.txt
#$ -o bash_output.txt
#$ -m ae
#$ -M nstorten@mpi-bremen.de


echo "job started:"
echo "job ID:$JOB_ID"
echo "date: $(date)"
echo "hostname: $(hostname)"


wdir=/scratch/tmp-$JOB_ID
mkdir $wdir

echo "------------and now python script"

python readQC_rawread.py $NSLOTS $wdir &> python_output.txt

echo "-----------that was the python script"

rm -rf $wdir
