args <- commandArgs(trailingOnly = TRUE)
    id   = as.character(args[1])
    files = args[2]
    vcffile = args[3]

vcf = read.table(vcffile) 
wigs = read.table(files) 
chrms = unique(matrix(unlist(strsplit( as.character(wigs[,1]), "\\.")), byr=T)[,1] ) 

 
list.vcf=list()
i = 1

for( chrm in chrms){
    f1 = vcf[,1]==chrm
    test1 = vcf[f1,]
    filewig =  paste( chrm,id,"intersect.wig", sep=".")
    info = file.info(filewig)

   if( !file.exists(filewig) ) { i = i + 1 ; next;  }  
   if( info$size==0 ) { i = i + 1 ; next;  } 
        
	test = read.table( filewig)
	test2 = test[,1:6]
	test2[,2:6] = test2[,2:6] + test[,9:13]
	rm(test)
	colnames(test2) = c( "Pos", "A", "C", "G", "T", "N" )
	sums = rowSums(test2[,2:6])
	b = rowSums(test2[,2:5]>0)

	f =  b >= 2

	if( sum(f) == 0) { i = i + 1 ; next; }


	f.all =  f
	d = test2[f.all,2:5]/sums[f.all]
	d[d==0]=NA
	d = as.matrix(d)
	test3.mask = test2[ f.all ,1:5]
	# test3.mask[test3.mask < 10] = 0
	b3.mask = rowSums(test3.mask[,2:5]>0) >=2 



	if( sum(f.all) == 0) { i = i + 1 ; next; }

	m = match(test3.mask[b3.mask,1], test1[,2])
	f.t = !is.na(m)
	f.v = m[f.t]

	if( sum(f.t) <= 1) { i = i + 1 ; next; }

	temp.vcf = cbind(d[b3.mask,][f.t,], test1[f.v,4:5])
	ref1 = as.character(temp.vcf[,5]) == "A"
	ref2 = as.character(temp.vcf[,5]) == "C"
	ref3 = as.character(temp.vcf[,5]) == "G"
	ref4 = as.character(temp.vcf[,5]) == "T"

	alt1 = as.character(temp.vcf[,6]) == "A"
	alt2 = as.character(temp.vcf[,6]) == "C"
	alt3 = as.character(temp.vcf[,6]) == "G"
	alt4 = as.character(temp.vcf[,6]) == "T"

	ref.only = c(temp.vcf[ref1,1], temp.vcf[ref2,2], temp.vcf[ref3,3], temp.vcf[ref4,4])
	alt.only = c(temp.vcf[alt1,1], temp.vcf[alt2,2], temp.vcf[alt3,3], temp.vcf[alt4,4])
	temp2.vcf = cbind(test3.mask[b3.mask,2:5][f.t,], test1[f.v,4:5])
	temp3.vcf = temp2.vcf[,1:2]
	temp3.vcf = temp3.vcf* 0

	temp3.vcf[ref1,1] = temp2.vcf[ref1,1]
	temp3.vcf[ref2,1] = temp2.vcf[ref2,2]
	temp3.vcf[ref3,1] = temp2.vcf[ref3,3]
	temp3.vcf[ref4,1] = temp2.vcf[ref4,4]

	temp3.vcf[alt1,2] = temp2.vcf[alt1,1]
	temp3.vcf[alt2,2] = temp2.vcf[alt2,2]
	temp3.vcf[alt3,2] = temp2.vcf[alt3,3]
	temp3.vcf[alt4,2] = temp2.vcf[alt4,4]

	x = c(alt.only, ref.only)
        list.vcf[[i]] = cbind(test3.mask[b3.mask,1][f.t], temp3.vcf, temp2.vcf)
        i = i + 1 
}

#list.vcf[[i-1]] = ""
#names(list.vcf) = chrms 
 skip = sapply( 1:length(list.vcf), function(i) length(list.vcf[[i]] ) )  ==0

#skip[i-1] = TRUE 
save(list.vcf, chrms, skip, file=paste0(id,".counts.vcf.Rdata"))

