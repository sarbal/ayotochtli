# Finding the X chromosome in the armadillo scaffolds 

### Mouse homologs 
load("homologs.Rdata")
load("gene_annotations_vM15.Rdata")
attr.mouse = attr
load("gene_annotations_v0.95_mod.Rdata")


chrX.genes = attr.mouse$name[attr.mouse$chr == "chrX"] 
m = match( genes.id2ids[,1], chrX.genes ) 
f.gi = !is.na(m)
f.xcg = m[f.gi]

m = match(genes.id2ids[f.gi,3], attr$name) 
f.x = !is.na(m)
f.ax = m[f.x]
freqX = count(attr$scaffold[f.ax])

Xscaffolds.mouse = as.character(freqX[,1])



### Human homologs
chain = read.table("chrX.chainHg38.txt")
scaff.sizes = unique(chain[,3:4]  )

library(rtracklayer)
# human 
chainObject <- import.chain("hg38ToDasNov3.over.chain")
grObject <- GRanges(seqnames=c("chrX"), ranges=IRanges(start=2, end=156040895))
results <- as.data.frame(liftOver(grObject, chainObject))
hist(results[,6], xlab="Length of alignment", main="hg38")

dim(results)
scaff.sums = (tapply( results[,6], results[,3], sum) )
scaff.range1 = (tapply( results[,4], results[,3], min) )
scaff.range2 = (tapply( results[,5], results[,3], max) )
scaff.range = scaff.range2 - scaff.range1

m = match(scaff.sizes[,1] , names(scaff.sums) )
f.s = !is.na(m)
f.ss = m[f.s]

prop = scaff.range[f.ss]/scaff.sizes[f.s,2]

plot(  scaff.sizes[f.s,2] , prop, pch=19 , xlab="Scaffold length (bp)", ylab="Fraction aligned to chrX")
which( prop > 0.8 & scaff.sizes[f.s,2] > 5e6)



m = match(scaffoldsX.sub, names(prop) )
f.sx = !is.na(m)
f.sr = m[f.sx]

plot(  log10(scaff.sizes[f.s,2]) , prop, pch=19 , xlab="Scaffold length (log10 bp)", ylab="Fraction aligned to chrX")
points(  log10(scaff.sizes[f.s,2][f.sr]) , prop[f.sr], pch=19 , col=2) 


zones=matrix(c(2,0,1,3), ncol=2, byrow=TRUE)
layout(zones, widths=c(4/5,1/5), heights=c(1/5,4/5))
x = log10(scaff.sizes[f.s,2])
y = prop  
xlab="Scaffold length (log10 bp)"
ylab= "Fraction aligned to chrX"
xhist0 = hist(x, plot=FALSE)
yhist0 = hist( y, plot=FALSE)
xhist1 = hist(x[f.sr], plot=FALSE, breaks=xhist0$breaks)
yhist1 = hist( y[f.sr], plot=FALSE, breaks=yhist0$breaks)
yhist1$counts = yhist1$counts/sum(yhist1$counts)
yhist0$counts = yhist0$counts/sum(yhist0$counts)
xhist1$counts = xhist1$counts/sum(xhist1$counts)
xhist0$counts = xhist0$counts/sum(xhist0$counts)

top = max(c(xhist1$counts, xhist0$counts, yhist0$counts, yhist1$counts))
par(mar=c(3,3,1,1))
plot(x,y, pch=19, main="", col=makeTransparent(1,50))
points( x[f.sr] , y[f.sr], pch=19 , col=2, cex=1.2) 

par(mar=c(0,3,1,1))
barplot(xhist0$counts, axes=FALSE, ylim=c(0, top), space=0)
axis(2)
barplot(xhist1$counts, axes=FALSE, ylim=c(0, top), space=0, add=T, col=makeTransparent(2)) 

par(mar=c(3,0,1,1))
barplot(yhist0$counts, axes=FALSE, xlim=c(0, top), space=0, horiz=TRUE)
axis(1)
barplot(yhist1$counts, axes=FALSE, xlim=c(0, top), space=0, horiz=TRUE, add=T, col=makeTransparent(2))

par(oma=c(3,3,0,0))
mtext(xlab, side=1, line=1, outer=TRUE, adj=0,
      at=.5 * (mean(x) - min(x))/(max(x)-min(x)))
mtext(ylab, side=2, line=1, outer=TRUE, adj=0,
      at=(.8 * (mean(y) - min(y))/(max(y) - min(y))))



# redo with alignment scaffolds 
frac = 0.9 
scaffoldsX.prop = names(prop[which(prop>= frac)])
scaffoldsX.prop = scaffoldsX.prop[-length(scaffoldsX.prop)] # remove chrM

frac = 0.9 
scaffoldsX.prop = names(prop[f.sr][which(prop[f.sr]>= frac)])

frac = 0.9 
scaffoldsX.prop2 = unique(c(names(which(prop > 0.99 & log10(scaff.sizes[f.s,2]) > 5 )) , names(prop[f.sr][which(prop[f.sr]>= frac)]))) 

scaffoldsX.prop3 =  c(scaffoldsX.prop2 , "JH573670")

frac = 0.9 
scaffoldsX.prop = unique(c(names(which(prop > 0.95 & log10(scaff.sizes[f.s,2]) > 5 )) , names(prop[f.sr][which(prop[f.sr]>= frac)]))) 


