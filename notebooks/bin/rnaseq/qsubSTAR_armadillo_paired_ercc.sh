dirFASTQ=$1
filesIn=`cd $dirFASTQ; ls *_1.fastq`
filesN=`echo $filesIn | wc -w`
nM=999

optEnc="--outFilterType BySJout   --outFilterMultimapNmax 20   --outFilterMismatchNoverReadLmax 0.04   --outFilterMatchNminOverLread 0 \
 --alignIntronMin 20   --alignIntronMax 1000000   --alignMatesGapMax 1000000   --alignSJoverhangMin 8   --alignSJDBoverhangMin 1   --sjdbScore 1 \
 --outFilterScoreMinOverLread 0.66  --outFilterMismatchNmax $nM  --outSAMtype BAM Unsorted  --quantMode GeneCounts  \
 --twopassMode Basic  --twopass1readsN -1 "

qsub -j y -v PATH -pe threads 10 -l tmp_free=200G,m_mem_free=20G -cwd -t 1-$filesN \
subSTAR_armadillo_paired_ercc.sh \
$dirFASTQ \"$filesIn\" ./ \
\" $optEnc \"
