#!/bin/bash
# run job in the current working directory where qsub is executed from
#$ -cwd
#  specify that the job requires 20GB of memory
#$ -l m_mem_free=20G
#$ -v OMP_NUM_THREADS

export OMP_NUM_THREADS=16
file=$1
bedfile=$2

perl ~/split_bed_wig_fast_updated.pl $file $bedfile
