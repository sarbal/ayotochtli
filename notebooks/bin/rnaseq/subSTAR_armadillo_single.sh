STAR=/sonas-hs/gillis/hpc/home/sballouz/STAR/STAR2
STAR=/sonas-hs/gillis/hpc/home/sballouz/STAR/STAR_2.5.2b
genome=/sonas-hs/gillis/hpc/home/sballouz/sballouz/armadillo/ensembl/Dasnov3.0.87/

inDir=$1
inRead1=( $2 )
checkDir=$1
STARoptions=$4

echo $checkDir
echo $STARoptions
echo $3

echo Number of files ${#inRead1[*]}
outPrefix=`pwd`

inI=$((SGE_TASK_ID-1))
read1=${inRead1[$inI]}

ls -l $inDir/$read1

readSuffix=${read1##*\.}

outDir=$inDir
echo $read1
echo $readSuffix $outDir $inDir

#mkdir $outDir
#cd $outDir
ls -lh $checkDir/$outDir/Log.final.out

if [ -s "$checkDir/$outDir/Log.final.out" ]
then
   echo "Run was completed"
   exit 0
fi

if [ "$readSuffix" == "fastq" ]
then
   readComm=-
   read1=$inDir$read1
fi


echo Read1: $read1

$STAR --genomeDir $genome --readFilesIn $read1 --readFilesCommand $readComm --runThreadN 5 $STARoptions

