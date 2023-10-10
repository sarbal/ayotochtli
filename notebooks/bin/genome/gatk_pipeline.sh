#!/bin/bash
# run job in the current working directory where qsub is executed from
#$ -cwd
#  specify that the job requires 20GB of memory
#$ -l m_mem_free=20G
#$ -v OMP_NUM_THREADS

export OMP_NUM_THREADS=16

samtools='/opt/hpc/bin/samtools'
in=$1
file=$2
genome=$3
picard='GenomeAnalysisTK/picard.jar'
GATK='GenomeAnalysisTK/GenomeAnalysisTK.jar'
igvtools='IGVTools/igvtools.jar'
# only run once
samtools faidx $genome
java -jar picard.jar CreateSequenceDictionary R=$genome O=$genome.dict

echo "Preparing BAM file for SNP calling"
if [ ! -e "$in/$file.Aligned.sorted.out.bai" ]
then
        $samtools index $in/$file.Aligned.sorted.bam
fi

$samtools view -b -q 10 $in/$file.Aligned.sorted.bam > $in/$file.Aligned.sorted.filt.bam

echo "Adding read groups to bam file"
java -jar $picard AddOrReplaceReadGroups \
    I=$in/$file.Aligned.sorted.filt.bam \
    O=$in/$file.Aligned.sorted.rg.bam \
    SO=coordinate RGID=id RGLB=library RGPL=platform RGPU=machine RGSM=sample

echo "Marking duplicates"
java -jar $picard MarkDuplicates \
    I=$in/$file.Aligned.sorted.rg.bam   \
    O=$in/$file.Aligned.sorted.dedupped.bam  \
    CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=output.metrics

echo "Splitting and trimming"
java -jar $GATK \
   -T SplitNCigarReads         \
   -R $genome \
   -I $in/$file.Aligned.sorted.dedupped.bam  \
   -o $in/$file.Aligned.sorted.split.bam  \
   -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS

 echo "Counting SNPs"
 java -jar $igvtools count -z 0 -w 1 --bases --strands read   \
   $in/$file.Aligned.sorted.split.bam   \
   $in/$file.filtered.wig   \
   $genome
