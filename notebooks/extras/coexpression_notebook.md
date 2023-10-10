
# Coexpression analysis
## Make armadillo aggregate network
```{r, eval=FALSE}

exprs = X.cpm.all
genes =rownames( X.cpm.all)
keep = rowSums(exprs) > 0


for( i in 1:n_quads){
  
   for (k in 1:n_times) {
     
      j=(k-1)*n_q;
      f.c = which(pData$Quad==i)[(1:n_q)+j]

      sub = exprs[keep,f.c]
      net = build_coexp_network(sub, genes[keep])
      med = median(net, na.rm=T)
      net[is.na(net)] = med
      
      # Save individual time point and quad ranked networks 
      save(net, file=paste("coexp.quad",i,"t", k, "Rdata", sep=".") )
      
      if(i==1 & k ==1 ) {
        agg = net
        } else {
          agg = net + agg
        }
      
   }
}

save(agg, file="coexp.agg.Rdata")

agg.rank =  matrix(rank(agg, na.last = "keep", ties.method = "average"), nrow=dim(agg)[1], ncol=dim(agg)[2] )
rownames(agg.rank) = rownames(agg)
colnames(agg.rank) = rownames(agg)
agg.rank = agg.rank/max(agg.rank, na.rm=T)
save(agg.rank, file="coexp.agg.rank.Rdata")


```


## Make quad aggregate networks
```{r, eval=FALSE}

for( i in 1:n_quads){
  print(i) 
  sub = sub * 0 
  agg = net * 0 
  for (k in 1:3) {
     print(k) 
      load(file=paste("coexp.quad",i,"t", k, "Rdata", sep=".") )
      agg = rowSums( cbind(agg, net), na.rm=T) 
  }
  print("agg")
  agg = rank(agg, na.last = "keep", ties.method = "average")
  agg = agg/max(agg, na.rm=T)
  save(agg, file=paste("coexp.agg.quad",i, "Rdata", sep=".") )
}   

```



## Make leave-one-quad-out aggregate networks and evaluate
```{r, eval=FALSE}
bottom = row(net) > col(net)

for( j in 1:n_quads){
  print(j) 
  sub = sub * 0 
  agg = sub * 0 
  loocv.quad = (1:n_quads)[-j]
  for( i in loocv.quad){
    print(i) 
     for (k in 1:n_times) {
       print(k) 
        load(file=paste("coexp.quad",i,"t", k, "Rdata", sep=".") )
        agg = rowSums( cbind(agg, net), na.rm=T) 
     }
  }
  print("agg")
  agg = rank(agg, na.last = "keep", ties.method = "average")
  agg = agg/max(agg, na.rm=T)
  save(agg, bottom,genes, f.zz, file=paste("coexp.quad.loo",j,"Rdata", sep=".") )  
}  
  

for( j in 1:n_quads){
  print(j) 
  sub = sub * 0 
  load(file=paste("coexp.quad.loo",j,"Rdata", sep=".") )  

  sub[bottom] = agg 
  sub = sub + t(sub)
  diag(sub) =  1 
  print("aurocs") 
  aurocs = run_GBA(sub, perfect.pred.clust, min=0)
  aurocs[[2]] = ""
  save(aurocs, file=paste("rocs.perfect.pred.clust.coexp.quad.loo",j,"Rdata", sep=".") )
}
```


## Make quad aggregate networks and evaluate
```{r, eval=FALSE}

for( i in 1:n_quads){
  print(i) 
  sub = sub * 0 
  agg = net * 0 
  for (k in 1:n_times) {
     print(k) 
      load(file=paste("coexp.quad",i,"t", k, "Rdata", sep=".") )
      agg = rowSums( cbind(agg, net), na.rm=T) 
  }
  print("agg")
  agg = rank(agg, na.last = "keep", ties.method = "average")
  agg = agg/max(agg, na.rm=T)
  save(agg, file=paste("coexp.agg.quad",i, "Rdata", sep=".") )
}   

for( i in 1:n_quads){
  print(i) 
  load(file=paste("coexp.agg.quad",i,"Rdata", sep=".") )
  sub = sub * 0   
  sub[bottom] = agg 
  sub = sub + t(sub)
  diag(sub) =  1 
  print("aurocs") 
  aurocs = run_GBA(sub, perfect.pred.clust, min=0)
  aurocs[[2]] = ""
  save(aurocs, file=paste("rocs.perfect.pred.clust.coexp.quad",i,"Rdata", sep=".") )
}

```

```{r}
aurocs.loo = list() 
summary.aurocs.clust.loo = list()

for(j in 1:n_quads){
load(file=paste("coexp/rocs.perfect.pred.clust.coexp.quad.loo",j,"Rdata", sep=".") ) 
  aurocs.loo[[j]] = aurocs[[1]][1:6,]
  summary.aurocs.clust.loo[[j]] = aurocs[[1]][-(1:6),]
}

aurocs.circ = list() 
summary.aurocs.clust.circ = list()

for(j in 1:n_quads){
load(file=paste("coexp/rocs.perfect.pred.clust.coexp.quad",j,"Rdata", sep=".") ) 
  aurocs.circ[[j]] = aurocs[[1]][1:6,]
  summary.aurocs.clust.circ[[j]] = aurocs[[1]][-(1:6),]
}

 

summary.aurocs.clust.circ.mat = summary.aurocs.clust.circ[[j]][,1:3]  * 0 
summary.aurocs.clust.loo.mat = summary.aurocs.clust.loo[[j]][,1:3]   * 0 
quadids = matrix(unlist(strsplit(rownames( summary.aurocs.clust.circ.mat  ), " " )  )  , ncol=2, by=T)[,1]

for(j in 1:n_quads){
  id = quads[j]
  summary.aurocs.clust.circ.mat[quadids==id,] = summary.aurocs.clust.circ[[j]][quadids==id, ]
  summary.aurocs.clust.loo.mat[quadids==id,] = summary.aurocs.clust.loo[[j]][quadids==id, ]
} 

filt  = (colSums(perfect.pred.clust)[-c(1:6)] > 1 )

aurocs.loo.list = lapply(1:n_quads, function(j) summary.aurocs.clust.loo.mat[quadids==quads[j] & filt,1] )
aurocs.dens = lapply(n_quads:1, function(i) density(aurocs.loo.list[[i]],   from=0, to=1)  )   

o = order(sapply(1:n_quads, function(i) mean( aurocs.loo.list[[i]] )), decreasing = F) 
aurocs.dens.sums = sapply(o, function(i) aurocs.dens[[i]]$y )
plot(aurocs.dens[[1]]$x, aurocs.dens[[1]]$y, xlim=c(0,1), ylim=c(0,10), col=0, xlab="AUROCs", ylab="Density") 
for(i in n_quads:2){ 
  polygon( c(0,aurocs.dens[[1]]$x,1), c(0,rowSums(aurocs.dens.sums[,1:i]),0) , col=(candy_colors[o])[i], border=NA)
} 
polygon( c(0,aurocs.dens[[1]]$x,1), c(0,aurocs.dens.sums[,1],0) , col=(candy_colors[o])[1], border=NA)
abline(v= mean( unlist(aurocs.loo.list[[i]] )), col=1, lwd=2)

```





## GO for armadillos
```{r}
load("GO.human.Rdata")
gosums = colSums(GO.human.nonIEA )
annot =  GO.human.nonIEA 

load("coexp/coexp.aggr.Rdata")
library(EGAD)

# Only use homologs
f.n = attr.arm$name[f.zz] != ""
network.arm =  sub[f.n,f.n]
bottom = row(network.arm) > col(network.arm)
sub = network.arm[bottom]
agg = rank(sub, na.last = "keep", ties.method = "average")
agg = agg/max(agg, na.rm=T)
sub = network.arm * 0
sub[bottom] = agg
sub = sub + t(sub)
diag(sub) = 1
rownames(sub) = attr.arm$name[f.zz]
colnames(sub) = attr.arm$name[f.zz]

aurocs.go1 = run_GBA( sub, annot,min=0 ,max=100 )
aurocs.go1[[2]] <- NULL 
save(aurocs.go1, file="aurocs.go1.arm.clust.Rdata")

```




## Human co-expression evaluation
```{r, eval=FALSE} 

load("gene_annotations_v29.Rdata")
attr.human = attr
load("gene_annotations_v0.95_mod.Rdata")
attr.arm = attr
### this file is too large, need to upload somewhere else
netfile="/sonas-hs/gillis/hpc/home/sballouz/sballouz/human/RNAseq/outliers/networks/blood.rerank.Rdata"
label="blood.rerank"
    nettype= label
load(netfile)

    network = diag(length(genes.t))
    bottom = row(network) > col(network)
    colnames(network) = genes.t
    rownames(network) = genes.t
    network[bottom] = temp
    network = network + t(network)
    diag(network) = 1


m = match( rownames(network), attr.human$entrezID )
f.n = !is.na(m)
f.ah = m[f.n]

network = network[f.n,f.n]
rownames(network) = attr.human$name[f.ah]
colnames(network) = attr.human$name[f.ah]


m2 = match( attr.human$name[f.ah], attr.arm$name[f.zz] ) 
f.n2 = !is.na(m2) 
f.aa2 = m2[f.n2]

network.arm =  network[f.n2,f.n2]

perfect.pred.human.clust = perfect.pred.clust[f.zz,][f.aa2,]
aurocs = run_GBA( network.arm, perfect.pred.human.clust,min=0 )




bottom = row(network.arm) > col(network.arm)
sub = network.arm[bottom]
agg = rank(sub, na.last = "keep", ties.method = "average")
agg = agg/max(agg, na.rm=T)
sub = network.arm * 0
sub[bottom] = agg 
sub = sub + t(sub) 
diag(sub) = 1 

aurocs = run_GBA(sub, perfect.pred.human.clust,min=0 )
aurocs[[2]] <- NULL
save(aurocs, perfect.pred.human.clust, file="perfect.pred.clust.human.bloodagg.reranked.Rdata")


annot = GO.human.nonIEA[f.g,]
rownames(annot) = attr.human$name[f.a] 


aurocs.go1 = run_GBA( sub, annot,min=0 ,max=100 )
aurocs.go1[[2]] <- NULL 

```

# Cross-species co-expression 
```{r}
hb = hist( unlist(aurocs.loo.list), breaks=10, freq=F) 
aurocs.hist = sapply(5:1, function(i) hist(aurocs.loo.list[[i]], breaks=hb$breaks, freq=F)$density )   
barplot(t(aurocs.hist/10), col=rev(candy_colors[1:5] ) , ylim=c(0,1.5), border=NA, space=0, names=hb$breaks[-1]) 



load("../coexp/perfect.pred.clust.human.bloodagg.reranked.Rdata")
filt = (colSums(perfect.pred.human.clust)[colSums(perfect.pred.human.clust) > 0] > 1  )
aurocs.loo.list$human = aurocs[[1]][filt,1]


aurocs.hum = sapply(5:1, function(i) hist(aurocs.loo.list$human , breaks=hb$breaks, freq=F)$density )   


load("../coexp/null_modules.blood.Rdata")
hist(aurocs.null, breaks=hb$breaks, col="lightgrey", border=NA, freq=F, main="", xlab="AUROCs")
abline(v=mean(aurocs.loo.list$human), lwd=2, col=2 ) 
abline(v=mean(aurocs.null), lty=2, lwd=2 ) 


load("")
bottom = row(sub.arm) > col(sub.arm)
armrank = rank(array(sub.arm[bottom]))
humanrank = rank(array(sub.blood[bottom]))

#pdf("temp.pdf")
#EGAD::conv_smoother( armrank, humanrank,  xlab="arm", ylab="human", 1000)
#dev.off() 


```


```{r}
load("coexp/aurocs.go1.blood.clust.Rdata") 
aurocs.go1.blood = aurocs.go1 
load("coexp/aurocs.go1.arm.clust.Rdata") 
aurocs.go1.arm = aurocs.go1 

m = match(rownames(aurocs.go1.blood[[1]]), rownames(aurocs.go1.arm[[1]] ) )
f.h = !is.na(m) 
f.a = m[f.h]

load("coexp/annot.agg.Rdata") 
mn = match( names(aurocs.go1.blood[[1]][f.h,1]), names(gosums )) 


main="GO"
xlab="GO AUROCs - human"
ylab="GO AUROCs - armadillo"


fgo = gosums[mn] >= 20 & gosums[mn] <= 1000

filt1_sig = rowSums((results.mat[filt1,4:8]) >= 4 ) > 0 
temp = GO.voc[rownames(results.mat[filt1,][filt1_sig,]),]


x =  aurocs.go1.blood[[1]][f.h,1][fgo]
y =  aurocs.go1.arm[[1]][f.a,1][fgo] 
sub = names(which(filt1_sig)) 
modcols = 0*(1:length(sub) )
for(i in 1:n_quads){ 
  modcols[ results.mat[filt1,i+3][filt1_sig] == 4  ]  = candy_colors[i]
} 

names(modcols) = names(which(filt1_sig))


plot(x,y, pch=19, col= "lightgrey", xlab=xlab, ylab=ylab) 
abline(0,1, lwd=2)
points(x[sub], y[sub], col=modcols, pch=19)
abline(v=mean(x[sub]), col=2, lwd=2 )
abline(h=mean(y[sub]), col=2, lwd=2 )
abline(v=mean(x), col="lightgrey", lty=2 )
abline(h=mean(y), col="lightgrey", lty=2 )

```




# Identity features: MsigDB
```{r}

load(file="msignrunlater.Rdata" )

xlab="MsigDB set size"
nr = 1:dim(genesetsgo)[2]
nr2 = colSums(genesetsgo[f.zz,]) 

pred.scores.list = list() 
pval.scores.list = list() 
pred.indscores.list = list()
setsizes.scores.list  = list() 
genesets.scores.list = list() 

i = 1 
for(k in which(nr2>5) ){ 
  pred.temp = matrix(0, ncol=3, nrow=5)
  pval.temp = matrix(0, ncol=3, nrow=5)
  pred.ind.temp = matrix(0, ncol=3*4, nrow=5)
  setsizes.temp = matrix(0, ncol=3, nrow=5) 
  genesets.temp = list() 
  
  for(j in 1:5){
    inds = pData$ID[c(((1:4)+4*(j-1)), ((1:4)+4*(j-1)) + 20, ((1:4)+4*(j-1)) + 40)]
    ki = ((1:4)+4*(j-1))
    kj = k*100
    
    q1 = genesetsgo[f.zz,k] == 1  
    q2 = q1
    q3 = q1
    
    t1 = prediction_gene_scores( X.tr1[f.zz,ki], X.tt1[f.zz,ki], q1) 
    t2 = prediction_gene_scores( X.tr2[f.zz,ki], X.tt2[f.zz,ki], q2) 
    t3 = prediction_gene_scores( X.tr3[f.zz,ki], X.tt3[f.zz,ki], q3) 
    
    pred.temp[j,] = c(t1[[2]],t2[[2]],t3[[2]])
    pred.ind.temp[j,match(c(t1[[3]],t2[[3]],t3[[3]]), inds)] = 1
    pval.temp[j,] = pval.raw[match(pred.temp[j,], pval.raw[,1]),2]
    
    setsizes.temp[j,] = c(length(q1), length(q2), length(q3)) 
    genesets.temp[[j]] = list((q1), (q2), (q3))   
  }
  pred.scores.list[[i]] = pred.temp
  pval.scores.list[[i]] = pval.temp
  pred.indscores.list[[i]] = pred.ind.temp
  setsizes.scores.list[[i]] = setsizes.temp
  genesets.scores.list[[i]] = genesets.temp
  
  i = i + 1
}



cd = sapply(nr, function(i) mean(pred.scores.list[[i]]))
bc = sapply(nr, function(i) colMeans(pred.scores.list[[i]]))
ab = sapply(nr, function(i) rowMeans(pred.scores.list[[i]]))
de = sapply(nr, function(i) (pred.scores.list[[i]]))
de.pval = sapply(nr, function(i) (pval.scores.list[[i]]))
bc.pval = matrix(pval.csum[match(bc, pval.csum[,1]),2], nrow=3)
ab.pval = matrix(pval.rsum[match(ab, pval.rsum[,1]),2], nrow=5)
cd.pval = pval.sum[match(cd, pval.sum[,1]),2]



save(
  genesetsgo, pred.scores.list,
  pval.scores.list,
  pred.indscores.list,
  setsizes.scores.list, 
  genesets.scores.list,nr, nr2, 
  file="task1_results_msig.Rdata")



#load("task1_results_msig.Rdata")
nr22 = nr2[nr2>5]
cd = sapply(1:length(nr22), function(i) mean(pred.scores.list[[i]]))
bc = sapply(1:length(nr22), function(i) colMeans(pred.scores.list[[i]]))
ab = sapply(1:length(nr22), function(i) rowMeans(pred.scores.list[[i]]))
de = sapply(1:length(nr22), function(i) (pred.scores.list[[i]]))
de.pval = sapply(1:length(nr22), function(i) (pval.scores.list[[i]]))
bc.pval = matrix(pval.csum[match(bc, pval.csum[,1]),2], nrow=3)
ab.pval = matrix(pval.rsum[match(ab, pval.rsum[,1]),2], nrow=5)
cd.pval = pval.sum[match(cd, pval.sum[,1]),2]
results.mat = cbind(nr22, cd,cd.pval, t(ab),t(ab.pval),  t(bc), t(bc.pval),  t(de), t(de.pval))
save(results.mat,file="task1_results_msig.mat.Rdata")


```






```{r} 
load(file="task1_results_msig.mat.Rdata")
filtk = grep("KEGG", rownames(results.mat) )
filtr = grep("REACTOME", rownames(results.mat) )
filtbc = grep("BIOCARTA", rownames(results.mat) )
filthlm = grep("HALLMARK", rownames(results.mat) )

filt_sig = rowSums((results.mat[filtk,4:8]) >= 4 ) > 0 
filt_sig = rowSums((results.mat[filthlm,4:8]) >= 3 ) > 0 
filt_sig = rowSums((results.mat[filtr,4:8]) >= 4 ) > 0 
filt_sig = rowSums((results.mat[filtbc,4:8]) >= 4 ) > 0 


heatmap.3( (results.mat[filtk,4:8][filt_sig ,]),Colv=F, col=magma(13), ColSideCol=candy_colors[1:5], labCol=quads)
heatmap.3( -log10(results.mat[filtk,9:13][ ,]),Colv=F, col=cols5, ColSideCol=candy_colors[1:5], labCol=quads)

heatmap.3( (results.mat[filthlm,4:8][  ,]),Colv=F, col=magma(13), ColSideCol=candy_colors[1:5], labCol=quads)
 
```




# Identity: tabula muris
```{r, eval=FALSE}
load("../updated/tab_mur_annot.Rdata")
m = match( attr$name, rownames(annot.tab) )
f.a = !is.na(m)
f.go = m[f.a]
genesetsgo = matrix(0, ncol=dim(annot.tab)[2], nrow=length(attr$name))
colnames(genesetsgo) = colnames(annot.tab)

for(i in 1:dim(annot.tab)[2]){ genesetsgo[f.a,i] = annot.tab[f.go,i]}
colSums(annot.tab)
colSums(genesetsgo)
genesetsgo = genesetsgo[,(colSums(genesetsgo[f.zz,]) >= 5 )]
genesetnamesgo = colnames(genesetsgo)


xlab="Tubula muris set size"
nr = 1:dim(genesetsgo)[2]
nr2 = colSums(genesetsgo[f.zz,]) 

pred.scores.list = list() 
pval.scores.list = list() 
pred.indscores.list = list()
setsizes.scores.list  = list() 
genesets.scores.list = list() 

i = 1 
for(k in which(nr2>5) ){ 
  pred.temp = matrix(0, ncol=3, nrow=5)
  pval.temp = matrix(0, ncol=3, nrow=5)
  pred.ind.temp = matrix(0, ncol=3*4, nrow=5)
  setsizes.temp = matrix(0, ncol=3, nrow=5) 
  genesets.temp = list() 
  
  for(j in 1:5){
    inds = pData$ID[c(((1:4)+4*(j-1)), ((1:4)+4*(j-1)) + 20, ((1:4)+4*(j-1)) + 40)]
    ki = ((1:4)+4*(j-1))
    kj = k*100
    
    q1 = genesetsgo[f.zz,k] == 1  
    q2 = q1
    q3 = q1
    
    t1 = prediction_gene_scores( X.tr1[f.zz,ki], X.tt1[f.zz,ki], q1) 
    t2 = prediction_gene_scores( X.tr2[f.zz,ki], X.tt2[f.zz,ki], q2) 
    t3 = prediction_gene_scores( X.tr3[f.zz,ki], X.tt3[f.zz,ki], q3) 
    
    pred.temp[j,] = c(t1[[2]],t2[[2]],t3[[2]])
    pred.ind.temp[j,match(c(t1[[3]],t2[[3]],t3[[3]]), inds)] = 1
    pval.temp[j,] = pval.raw[match(pred.temp[j,], pval.raw[,1]),2]
    
    setsizes.temp[j,] = c(length(q1), length(q2), length(q3)) 
    genesets.temp[[j]] = list((q1), (q2), (q3))   
  }
  pred.scores.list[[i]] = pred.temp
  pval.scores.list[[i]] = pval.temp
  pred.indscores.list[[i]] = pred.ind.temp
  setsizes.scores.list[[i]] = setsizes.temp
  genesets.scores.list[[i]] = genesets.temp
  
  i = i + 1
}



cd = sapply(nr, function(i) mean(pred.scores.list[[i]]))
bc = sapply(nr, function(i) colMeans(pred.scores.list[[i]]))
ab = sapply(nr, function(i) rowMeans(pred.scores.list[[i]]))
de = sapply(nr, function(i) (pred.scores.list[[i]]))
de.pval = sapply(nr, function(i) (pval.scores.list[[i]]))
bc.pval = matrix(pval.csum[match(bc, pval.csum[,1]),2], nrow=3)
ab.pval = matrix(pval.rsum[match(ab, pval.rsum[,1]),2], nrow=5)
cd.pval = pval.sum[match(cd, pval.sum[,1]),2]



save(
  genesetsgo, pred.scores.list,
  pval.scores.list,
  pred.indscores.list,
  setsizes.scores.list, 
  genesets.scores.list,nr, nr2, 
  file="task1_results_tab_mur.Rdata")

results.mat = cbind(nr, cd,cd.pval, t(ab),t(ab.pval),  t(bc), t(bc.pval),  t(de), t(de.pval))

save(results.mat, file="task1_results_tab_mur.mat.Rdata") 

f.ms= meta.gs.mat$hutchins_annotation == "Blood mesoderm"
heatmap.3(aurocs.tab[f.ms,], col=magma(100), ColSideCol=candy_colors[pData$Quad])


m = match( rownames(aurocs.tab[f.ms,] ), genesetnamesgo )
ff = !is.na(m)
f.ms2 = m[ff]

cellypes.col = sample(rainbow(50))[celltypes[f.ms2]] 
rownames(results.mat ) = genesetnamesgo

 
m = match( rownames(aurocs.tab[f.ms,] ), genesetnamesgo )
ff = !is.na(m)
f.ms2 = m[ff]

cellypes.col = sample(rainbow(50))[celltypes[f.ms2]] 
rownames(results.mat ) = genesetnamesgo

heatmap.3( (results.mat[f.ms2,4:8][ ,]),Colv=F, col=magma(13), ColSideCol=candy_colors[1:5], labCol=quads, 
           RowSideCol=cellypes.col) 
heatmap.3( -log10(results.mat[f.ms2,9:13][ ,]),Colv=F, col=cols5, ColSideCol=candy_colors[1:5], labCol=quads,
           RowSideCol=cellypes.col) 
```





# Comparison to historical phenotypic data
```{r}
load("data_from_old_papers/table1_scutes_newman1913.Rdata")
o = order(scutes[,2] ) 
pdf("data_from_old_papers/historicaldata.pdf")
boxplot( t(scutes[o,3:6]) ~ scutes[o,2] ) 
hist((rowSD(scutes[o,3:6]))^2, breaks=20) 
abline( v = sd( unlist(scutes[,3:6]))^2, lwd=2, col=2 )
abline( v = mean(rowSD(scutes[o,3:6])^2), lwd=2, col=4 )

load("data_from_old_papers/storrs_1968.Rdata")
boxplot( X$Body.weight ~ X$Set  )

boxplot( X$Scutes ~ X$Set  )
dev.off() 

```



# Compositional analysis
```{r}
frac = read.table("CIBERSORTx_Job5_Adjusted.txt", header=T, sep="\t")
frac2 = t(frac[,2:7])
barplot(frac2[,pData$Quad==j], col=inferno(6)) 
colnames(frac2) = pData$ID
pdf("cibersort.pdf")
for(j in 1:5) {
   heatmap.3( frac2[,pData$Quad==j] , 
           col= inferno(100), 
           ColSideCol=tropical_colors[pData$Ind[pData$Quad==j] ] ) 

 heatmap.3( frac2[,pData$Quad==j][,order(pData$altID[which(pData$Quad==j)] )] , 
           col= inferno(100), 
           ColSideCol=tropical_colors[pData$Ind[pData$Quad==j] ][order(pData$altID[which(pData$Quad==j)] )] , Colv=F)
} 
dev.off() 

```


