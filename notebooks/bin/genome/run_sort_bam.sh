dir=`pwd`
files=$1
for file in `cat $files | cut -f1 `
do
	echo $file
	cd bams
	qsub ../sort_bam.sh $dir/bams/$file
	cd ../
done
