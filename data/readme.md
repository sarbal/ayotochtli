### Data
#### Perfect predictor genes 
Genes that were able to predict armadillos across time, per quad can be found in this file: 
- perfect.pred.arm.Rdata 
This contains a matrix (perfect.pred) of dimensions 33374 (genes) by 5 (quads).  

The human and mouse homologs of thees genes are in these files: 
- perfect.pred.human.Rdata
- perfect.pred.mouse.Rdata

Each contains these: 
perfect.pred.<species> -> matrix, genes by quad
perfect.pred.<species>.list -> list of genes for each of the five quads 
homol.all2 -> gene conversion ids between species 

### Code
#### Perfect predictors
```
load("ref.strand.test_train.Rdata")
hist(pred.all[f.zz,], breaks=c(-1:4), col=magma(5), main="", xlab="Gene scores", sub=dataset)
hist(predcor.all[f.zz,], col=viridis(20), main="Testing", sub=dataset, xlab="Correlations")
hist(cor.all[f.zz,], col=viridis(20), main="Training", sub=dataset, xlab="Correlations")

perfect.pred.list = lapply( c(1,4,7,10,13), function(i) which(cor.all[f.zz,i]==1 & predcor.all[f.zz,i]==1 ) )
perfect.pred = sapply( c(1,4,7,10,13), function(i) (cor.all[ ,i]==1 & predcor.all[ ,i]==1 ) )

```
