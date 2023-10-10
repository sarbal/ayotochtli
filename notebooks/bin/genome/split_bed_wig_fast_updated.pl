#! /usr/bin/perl

#@HASH{@BARCODES} = 1;
$id =  $ARGV[0];
$bedfile =  $ARGV[1];
$file = "../wigs/".$id.".filtered.wig";

# $end=10006840167;



#`grep chrom $file -n > temp.chrm`;
#`cut -f1 $bedfile  | sort | uniq > temp.chrm2`;
#`cut -f2 -d' ' temp.chrm | sed 's/chrom=//g' | sort | uniq >  temp.chrm3`;
#@CHRMS = `grep -Fw -f temp.chrm2 temp.chrm3`;
#@HASH{@CRHMS} = 1; 

open (IN, "< $file ");
while(<IN>){
   chomp $_ ;
  if ( m/^variableStep chrom=(.+) span=1/ ) {
    $prev = $chrm;  
    $chrm = $1;
    close OUT; 
    
   `grep ^$prev $bedfile | cut -f2 > $prev.$id.temp.B`;
   `grep -Fw -f $prev.$id.temp.B $prev.$id.file >  $prev.$id.intersect.wig`;
   `rm  $prev.$id.temp.B  $prev.$id.file `;
    
    open( OUT, "> $chrm.$id.file");
 
 } else { 
  print OUT "$_\n";
 } 
}
close IN;


