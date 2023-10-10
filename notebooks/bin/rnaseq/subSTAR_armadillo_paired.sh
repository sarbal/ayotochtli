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
read2=${read1/_1/_2}
#read2=${read1/R1/R2}
#read2=${read2/R1/R2}

echo $read2

ls -l $inDir/$read1
ls -l $inDir/$read2

readSuffix=${read1##*\.}

outDir=$read1
#outDir=${outDir/.fastq.filtered_pairs_R1.fastq/""}
#outDir=${outDir/fastq.filtered_pairs_R1.fastq/ercc}
echo $read1 $readSuffix $outDir

mkdir $outDir
cd $outDir
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
   read2=$inDir$read2
fi

if [ "$readSuffix" == "gz" ]
then
   readComm=zcat
   read1=$inDir$read1
   read2=$inDir$read2
fi

if [ "$readSuffix" == "tgz" ]
then
   echo `date` : decompressing tgz files
   mkdir $TMPDIR/R1
   mkdir $TMPDIR/R2
   tar xfz $inDir$read1 --directory=$TMPDIR/R1
   tar xfz $inDir$read2 --directory=$TMPDIR/R2
   read1=`find $TMPDIR/R1 -name "*.*" | sort | tr '\n' ','`
   read2=`find $TMPDIR/R2 -name "*.*" | sort | tr '\n' ','`
   readComm="cat"
   echo `date` : done
fi

echo Read1: $read1
echo Read2: $read2

$STAR --genomeDir $genome --readFilesIn $read1 $read2 --readFilesCommand $readComm --runThreadN 1 $STARoptions

