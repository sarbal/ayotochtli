#!/bin/bash
# run job in the current working directory where qsub is executed from
#$ -cwd
#  specify that the job requires 20GB of memory
#$ -l m_mem_free=20G
#$ -v OMP_NUM_THREADS

export OMP_NUM_THREADS=16

samtools='/opt/hpc/bin/samtools'
file=$1

$samtools sort -@ 10 $file.Aligned.out.bam $file.Aligned.sortedByCoord.out
