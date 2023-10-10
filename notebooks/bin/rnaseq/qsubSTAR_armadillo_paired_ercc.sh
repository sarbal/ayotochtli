dirFASTQ=$1
filesIn=`cd $dirFASTQ; ls *_1.fastq`
filesN=`echo $filesIn | wc -w`
nM=999

optEnc="--outFilterType BySJout  --outSAMtype BAM Unsorted  --quantMode GeneCounts  \
 --twopassMode Basic  --twopass1readsN -1 --outFilterScoreMinOverLread 0.5 --outFilterMatchNminOverLread 0.5 "

qsub -j y -v PATH -pe threads 10 -l tmp_free=200G,m_mem_free=10G -cwd -t 1-$filesN \
 subSTAR_armadillo_paired_ercc2.sh \
$dirFASTQ \"$filesIn\" ./ \
\" $optEnc \"
