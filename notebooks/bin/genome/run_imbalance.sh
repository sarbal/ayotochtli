#!/bin/bash
# run job in the current working directory where qsub is executed from
#$ -cwd
#  specify that the job requires 20GB of memory
#$ -l m_mem_free=20G
#$ -v OMP_NUM_THREADS

export OMP_NUM_THREADS=16
id=$1
file=$2
bedfile=$3
bin='/armadillo/genomes/rna_fastq/bin/'

R  --no-save --args $id $file $bedfile < $bin/parse_imbalance.r 
