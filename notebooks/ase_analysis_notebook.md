---
title: "Armadillo transcritional identity - ASE"
output: html_notebook
---
# ASE data 
```{r}
source("armadillo_helper.r")
source("load_ase_data.r")
```
 
 
```{r, eval=FALSE}

exprs.all.filt = lapply( 1:5, function(i) exprs.all[[i]][!is.na(exprs.all[[i]][,1]),]) 
ratios = list() 
for(j in 1:5){ 
  X.temp = (exprs.all.filt[[j]][,7:18]) 
  nj = dim(X.temp)[1]/2
  ratios[[j]] = cbind((X.temp[(1:nj),1:4] / (X.temp[(1:nj)+nj,1:4] +X.temp[(1:nj),1:4] ) ),
                      (X.temp[(1:nj),5:8] / (X.temp[(1:nj)+nj,5:8] +X.temp[(1:nj),5:8] ) ),
                      (X.temp[(1:nj),9:12] / (X.temp[(1:nj)+nj,9:12] +X.temp[(1:nj),9:12] ) ))
  
}
density.plot = lapply(1:5, function(i) density( ratios[[i]][!is.na(ratios[[i]])] )  )



```


```{r, eval = FALSE}


## Min/max powered SNPs
exprs.all.filt.min = list() 
exprs.all.filt.max = list() 
for(j in 1:5){ 
  a = rowSums(exprs.all.filt[[j]][,7:18])
  b = exprs.all.filt[[j]][,1]
  nj = length(b)/2 
  
  abi = tapply(a,b, which.max)
  abi = abi[!is.na(abi)]
  
  abj = tapply(a,b, which.min)
  abj = abj[!is.na(abj)]
  
  aji = tapply(a,b, length)/2
  aji = aji[!is.na(aji)]
  
  abi_alt = abi+aji
  abi_alt[(abi-aji)>0] = (abi-aji)[(abi-aji)>0] 
  abj_alt = abj+aji 
  abj_alt[(abj-aji)>0] = (abj-aji)[(abj-aji)>0]
  
  abi = cbind(names(abi), abi, abi_alt)
  abj = cbind(names(abj), abj, abj_alt)
  
  
  maxtest  = lapply(1:dim(abi)[1], function(i) exprs.all.filt[[j]][ (exprs.all.filt[[j]][,1]==abi[i,1]),][ as.numeric(abi[i,2:3]), ] )
  mintest  = lapply(1:dim(abj)[1], function(i) exprs.all.filt[[j]][ (exprs.all.filt[[j]][,1]==abj[i,1]),][ as.numeric(abj[i,2:3]), ] )
  
  exprs.all.filt.max[[j]] = do.call(rbind, maxtest)
  exprs.all.filt.min[[j]] = do.call(rbind, mintest)
  
  
}



ratios.max.genes =  list() 
ratios.max = list() 
for(j in 1:5){ 
  X.temp = (exprs.all.filt.max[[j]][,7:18]) 
  nj = dim(X.temp)[1]/2
  ref = 2*(1:nj) - 1 
  alt = 2*(1:nj)  
  ratios.max.genes[[j]] =  cbind(exprs.all.filt.max[[j]][ref,1:6], exprs.all.filt.max[[j]][alt,6])
  ratios.max[[j]] = cbind((X.temp[ref,1:4] / (X.temp[alt,1:4] + X.temp[ref,1:4] ) ),
                          (X.temp[ref,5:8] / (X.temp[alt,5:8] + X.temp[ref,5:8] ) ),
                          (X.temp[ref,9:12]/ (X.temp[alt,9:12]+ X.temp[ref,9:12] ) ))
}


ratios.min = list() 
ratios.min.genes = list() 
for(j in 1:5){ 
  X.temp = (exprs.all.filt.min[[j]][,7:18]) 
  nj = dim(X.temp)[1]/2
  ref = 2*(1:nj) - 1 
  alt = 2*(1:nj)  
  ratios.min.genes[[j]] =  cbind(exprs.all.filt.min[[j]][ref,1:6], exprs.all.filt.min[[j]][alt,6])
  ratios.min[[j]] = cbind((X.temp[ref,1:4] / (X.temp[alt,1:4] + X.temp[ref,1:4] ) ),
                          (X.temp[ref,5:8] / (X.temp[alt,5:8] + X.temp[ref,5:8] ) ),
                          (X.temp[ref,9:12]/ (X.temp[alt,9:12]+ X.temp[ref,9:12] ) ))
  
  
}

ratios  = list() 
ratios.genes = list() 
for(j in 1:5){ 
  X.temp = (exprs.all.filt[[j]][,7:18]) 
  nj = dim(X.temp)[1]/2
  ref = 2*(1:nj) - 1 
  alt = 2*(1:nj)  
  ratios.genes[[j]] =  cbind(exprs.all.filt[[j]][ref,1:6], exprs.all.filt[[j]][alt,6])
  ratios[[j]] = cbind((X.temp[ref,1:4] / (X.temp[alt,1:4] + X.temp[ref,1:4] ) ),
                      (X.temp[ref,5:8] / (X.temp[alt,5:8] + X.temp[ref,5:8] ) ),
                      (X.temp[ref,9:12]/ (X.temp[alt,9:12]+ X.temp[ref,9:12] ) ))
  
  
}




# ratios
X.tt1 = list() 
X.tt2 = list() 
X.tt3 = list() 
X.tr1 = list() 
X.tr2 = list() 
X.tr3 = list()


cor.t1t2 = list()
cor.t1t3 = list()
cor.t2t3 = list()


for(j in 1:5){
  
  X.temp = apply(ratios.max[[j]], 2, as.numeric) 
  colnames(X.temp) = colnames(exprs.all.filt.max[[j]][,7:18]) 
  X.rank = lapply(1:3, function(i)   t(apply(X.temp[,(1:4) + 4*(i-1) ],  1, rank )))
  X.rank2 = do.call(cbind, X.rank)
  
  X.temp = apply(exprs.all.filt.max[[j]][,7:18], 2, as.numeric) 
  colnames(X.temp) = colnames(exprs.all.filt.max[[j]][,7:18]) 
  
  # 
  N = dim(X.temp)[1]
  
  # Testing 
  X.tt1[[j]] = X.rank2[,(1:4)]
  X.tt2[[j]] = X.rank2[,(5:8)]
  X.tt3[[j]] = X.rank2[,(9:12)]
  
  # Training
  X.t1 = X.temp[,-(1:4)]
  X.t2 = X.temp[,-(5:8)]
  X.t3 = X.temp[,-(9:12)]
  
  
  X.tr1q =  t(sapply(1:N, function(i) (rowMeans(cbind((X.t1[i, 1:4]), (X.t1[i, (1:4)+4 ])))))) 
  X.tr2q =  t(sapply(1:N, function(i) (rowMeans(cbind((X.t2[i, 1:4]), (X.t2[i, (1:4)+4 ]))))) )
  X.tr3q =  t(sapply(1:N, function(i) (rowMeans(cbind((X.t3[i, 1:4]), (X.t3[i, (1:4)+4 ]))))))  
  
  nj = N/2
  ref = 2*(1:nj) - 1 
  alt = 2*(1:nj)  
  
  X.tr1r =   X.tr1q[ref,1:4] / (X.tr1q[alt,1:4] + X.tr1q[ref,1:4] ) 
  X.tr2r =   X.tr2q[ref,1:4] / (X.tr2q[alt,1:4] + X.tr2q[ref,1:4] ) 
  X.tr3r =   X.tr3q[ref,1:4] / (X.tr3q[alt,1:4] + X.tr3q[ref,1:4] ) 
  
  X.tr1[[j]] = t(apply(X.tr1r, 1, rank ))
  X.tr2[[j]] = t(apply(X.tr2r, 1, rank )) 
  X.tr3[[j]] = t(apply(X.tr3r, 1, rank )) 
  
  
}




for(j in 1:5){
  
  X.temp = apply(ratios.max[[j]], 2, as.numeric) 
  colnames(X.temp) = colnames(exprs.all.filt.max[[j]][,7:18]) 
  X.rank = lapply(1:3, function(i)   t(apply(X.temp[,(1:4) + 4*(i-1) ],  1, rank )))
  X.rank2 = do.call(cbind, X.rank)
  
  
  
  N = dim(X.temp)[1]
  # Correlations between t1 and t2
  cor.t1t2[[j]] =  sapply(1:N, function(i) cor( (X.temp[i,1:4]), (X.temp[i,5:8]), m="s") )  
  # Correlations between t1 and t3
  cor.t1t3[[j]] =  sapply(1:N, function(i) cor( (X.temp[i,1:4]), (X.temp[i, 9:12]), m="s") )  
  # Correlations between t2 and t3 
  cor.t2t3[[j]] =  sapply(1:N, function(i) cor( (X.temp[i,5:8]), (X.temp[i, 9:12]), m="s") )  
  
}


rm(mintest)
rm(maxtest)

save( 
  exprs.all.filt, ratios.genes, 
  exprs.all.filt.max, ratios.max, ratios.max.genes, 
  exprs.all.filt.min, ratios.min, ratios.min.genes, 
  X.tt1, X.tt2, X.tt3, X.tr1, X.tr2, X.tr3,cor.t1t2,
  cor.t1t3, cor.t2t3, file= "ase_ratios.test_train.Rdata"  )


cor.all.list  = lapply(1:n_quads, function(j) cbind( cor.t2t3[[j]], cor.t1t3[[j]], cor.t1t2[[j]])) 
perfect.pred.list = lapply(1:n_quads, function(i) (rowSums(cor.all.list[[i]] == 1, na.rm=T  ) == 3 )*1 )  
perfect.pred.ase = perfect.pred * 0 
rownames(perfect.pred.ase) = attr$ensemblID
for(i in 1:5){ 
  m = match( attr$ensemblID, ratios.max.genes[[i]]$ensemblID)
  f.a = !is.na(m)
  f.r = m[f.a]
  perfect.pred.ase[f.a,i] =  perfect.pred.list[[i]][f.r]
}

cor.ase = cor.all  * NA 
rownames(cor.ase) = attr$ensemblID
for(i in 1:15){ 
  j = round(((i)+1)/3) 
  k = i %% 3 
  if( k == 0) { k = 3 }
  m = match( attr$ensemblID, ratios.max.genes[[j]]$ensemblID)
  f.a = !is.na(m)
  f.r = m[f.a]
  cor.ase[f.a,i] = cor.all.list[[j]][,k]   
}


f.na = sapply(1:15, function(i) !is.na(cor.ase[,i]   ) ) 

aurocs.comp = sapply(1:15, function(i) auroc_analytic( rank(cor.ase[f.zz&f.na[,i],i]), 1*(cor.all[f.zz&f.na[,i],i]==1) ) )  
aurocs.comp = sapply(1:15, function(i) auroc_analytic( rank(cor.ase[f.zz&f.na[,i],i]), 1*(cor.all[f.zz&f.na[,i],i]>0) ) )  

plots.smooth = sapply(1:15, function(i) EGAD::conv_smoother( cor.ase[f.zz&f.na[,i],i], cor.all[f.zz&f.na[,i],i], 50, "", "") ) 

plots.box = sapply(1:15, function(i) plot( cor.ase[f.zz&f.na[,i],i] ~ 1*(cor.all[f.zz&f.na[,i],i] > 0)  ) ) 

```

# SNPs 
```{r}
load("snps.overlaps.Rdata")
names(x)  = quads 
library(UpSetR)
upset(fromList(x), order.by = "freq", sets=names(x), keep.order=T, sets.bar.color = candy_colors[1:5])
library(venn)
venn( x, zcol=candy_colors , cexil = 1, cexsn = 1.5, ellipse=T)

```

```{r}

chrminfo = read.table("chromInfo.txt.gz")
chrmadd = cumsum(round(chrminfo[,2]/1e3))
chrmadd = c(0, chrmadd[-length(chrmadd)])
chrminfo[,3] = chrmadd*1e3

possnps = lapply(1:5, function(i) matrix(unlist(strsplit(x[[i]], " ")) , ncol=2, by=T) )
m = lapply(1:5, function(i) match(possnps[[i]][,1], chrminfo[,1]) )
possnps2 = lapply(1:5, function(i) as.numeric(possnps[[i]][,2]) + as.numeric(chrminfo[m[[i]],3]) )


hb = hist( round(unlist(possnps2)/1e3), breaks=1000, freq=F)

snpdensity = sapply(1:5, function(i) hist( round(possnps2[[i]]/1e3), breaks=hb$breaks, freq=F)$density )
heatmap.3(t(snpdensity), col=cividis(100), Colv=F)


```

# Genes with variants 
```{r}
x = lapply(1:5, function(i) ratios.max.genes[[i]][,1] )   
names(x)  = quads 
library(UpSetR)
upset(fromList(x), order.by = "freq", sets=names(x), keep.order=T, sets.bar.color = candy_colors[1:5])
library(venn)
venn( x, zcol=candy_colors , cexil = 1, cexsn = 1.5, ellipse=T)

``` 
 
# Allelic ratios
```{r}
density.plot.max = lapply(1:5, function(i) density( ratios.max[[i]][!is.na(ratios.max[[i]])], bw=0.03)  )
density.plot.min = lapply(1:5, function(i) density( ratios.min[[i]][!is.na(ratios.min[[i]])], bw=0.03)  )
density.plot = lapply(1:5, function(i) density( ratios[[i]][!is.na(ratios[[i]])], bw=0.03)  )
pdf("U:/armadillo/updated/density.ase.plot.pdf")
plot( density.plot.max[[1]], col=0, xlim=c(0,1), ylim=c(0,4.4), main="", xlab="Allelic ratios")
sapply(1:5, function(i) lines( density.plot.max[[i]], lwd=3, lty=2, col = makeTransparent(candy_colors[i]) ) )
sapply(1:5, function(i) lines( density.plot.min[[i]], lwd=3, lty=3, col = makeTransparent(candy_colors[i]) ) )
sapply(1:5, function(i) lines( density.plot[[i]], lwd=3, col = makeTransparent(candy_colors[i]) ) )
dev.off() 

```



# XCI
## Ratios 
```{r}
scaffoldsX = scaffoldsX.prop
# load(file="skew.est.max.genes.Rdata" )

folded <- function(x) { apply( cbind(x,1-x), 1, max)} 
unfold <- function(x) { c(x,1-x) } 

mle_folded <- function(){ 
  mus = seq(0.5,1, by = 0.001)
  sigmas = seq(0,0.5, by = 0.01)
  
  mles = sapply(mus, function(mu)
    sapply(sigmas, function(sigma)
      -sum( log(dnorm( x,   mean = mu, sd = sigma )
                + dnorm( x,   mean = 1-mu, sd = sigma ) )  )))
  mles[!is.finite(mles)] = NA
  coefs = which(mles==min(mles, na.rm=T), arr.ind=T)
  return (list( mus[coefs[2]] , sigmas[coefs[1]]))
} 



 
# max only 
exprs.maybeX.prop= list() 
ratios.maybeX.prop = list() 

for(j in 1:5){ 
  m = match(  exprs.all.filt.max[[j]]$chrm, scaffoldsX )
  f.imx = !is.na(m)
  f.i = m[f.imx]
  
  exprs.maybeX.prop[[j]] = exprs.all.filt.max[[j]][f.imx,] 
  
  m = match(  ratios.max.genes[[j]]$chrm, scaffoldsX )
  f.imx = !is.na(m)
  f.i = m[f.imx]
  
  ratios.maybeX.prop[[j]] = ratios.max[[j]][f.imx,]
} 


  
combined.ratios.list = list() 
folded.x.list = list() 
fittfold.comb.list = list()
  for(j in 2:4){    
    folded.x.list[[j]] = list()
    fittfold.comb.list[[j]] = list()
    nj = dim(exprs.maybeX.prop[[j]])[1]/2
    ref_counts = exprs.maybeX.prop[[j]][2*(1:nj)-1,  6 + 1:4 ] 
                +  exprs.maybeX.prop[[j]][2*(1:nj)-1,  6  + 5:8]
                +  exprs.maybeX.prop[[j]][2*(1:nj)-1, 6 + 9:12 ]
    
    alt_counts = exprs.maybeX.prop[[j]][2*(1:nj),  6 + 1:4 ] 
    +  exprs.maybeX.prop[[j]][2*(1:nj),  6  + 5:8]
    +  exprs.maybeX.prop[[j]][2*(1:nj), 6 + 9:12 ]
    

     combined.ratios.list[[j]] = ref_counts/(ref_counts+alt_counts) 
     
     for( i in 1:4){ 
        filt1 =  ref_counts[,i] > 5 &  alt_counts[,i] > 5
        x = folded(combined.ratios.list[[j]][filt1,i])
        folded.x.list[[j]][[i]] = x 
        fittfold = mle_folded()
        fittfold.comb.list[[j]][[i]] = fittfold
    } 
  }
  



est_fold_means = sapply(2:4, function(j) sapply(1:4, function(i) fittfold.comb.list [[j]][[i]][[1]])) 


# save(fittfold.comb.list,est_fold_means, est_fold_sds, folded.x.list , exprs.maybeX.prop, combined.ratios.list, file="estimated_ratios.Rdata" )



 # save(pvals.list, fittfold.list, est_fold.list, exprs.maybeX.prop, exprs.maybeX.sub, exprs.maybeX, ratios.maybeX.prop,ratios.maybeX, ratios.maybeX.sub,   file="skew.est.max.genes.Rdata" )
  





hist( unfold(est_fold_means), xlab="Estimated skew",  main="", col=makeTransparent(magma(10)[3]), border=NA, xlim=c(0,1), ylim=c(0,10))
lines( sort(unfold(p)), dnorm(sort(unfold(p)), mean = 0.5 , sd = sqrt(0.25/2) ) ,   lwd=2, col=magma(10)[1] ) 
lines( sort(unfold(p)), dnorm(sort(unfold(p)), mean = 0.5 , sd = sqrt(0.25/4) ) ,   lwd=2, col=magma(10)[2] ) 
lines( sort(unfold(p)), dnorm(sort(unfold(p)), mean = 0.5 , sd = sqrt(0.25/8)) ,    lwd=2, col=magma(10)[3] ) 
lines( sort(unfold(p)), dnorm(sort(unfold(p)), mean = 0.5 , sd = sqrt(0.25/16) ) ,  lwd=2, col=magma(10)[4] ) 
lines( sort(unfold(p)), dnorm(sort(unfold(p)), mean = 0.5 , sd = sqrt(0.25/32) ) ,  lwd=2, col=magma(10)[5] ) 
lines( sort(unfold(p)), dnorm(sort(unfold(p)), mean = 0.5 , sd = sqrt(0.25/64) ) ,  lwd=2, col=magma(10)[6] ) 
lines( sort(unfold(p)), dnorm(sort(unfold(p)), mean = 0.5 , sd = sqrt(0.25/100) ) , lwd=2, col=magma(10)[7] )  

lines( sort(unfold(p)), dnorm(sort(unfold(p)), mean = 0.5 , sd = 0.1 ) , lwd=2, col=magma(10)[7], lty=2 )  
abline(v=est, col=magma(10)[7],lwd=2, lty=2) 


  
```
##  Cell estimates 
```{r}
foo <- function(i,j) { 
  return( 0.25/(sd(unfold(i[j]))^ 2))
} 

est_from = tail(sort(array(est_fold_means)), n=3)

library(boot)
myBootstrap3 <- boot( est_from, foo  , R=1000 )
bootCI3 = boot.ci(myBootstrap3, index=1)

N = 64 
pi = 0.5
var_est = (pi * (1-pi)) / (1:N)
plot(  (1:N), (var_est), pch=19, xlab="N cells", ylab="Var")
abline(h=var, col=2, lwd=2)
abline(v= (nn), col=2,lwd=2)
polygon(  ( c(bootCI3$percent[4],bootCI3$percent[5],bootCI3$percent[5],bootCI3$percent[4])), 
          c(0,0,1,1),
          col=makeTransparent(magma(5)[4]),
          border=NA)

```



# ASE identity tests
## X scaffolds 
```{r}
# load("Xscaffolds.Rdata")
set.seed(5)
scaffoldsX = scaffoldsX.prop3 
scaffoldsX = scaffoldsX.sub

get_ase_score <- function(){
  pred.temp = matrix(0, ncol=3, nrow=5)
  pval.temp = matrix(0, ncol=3, nrow=5)
  pred.ind.temp = matrix(0, ncol=3*4, nrow=5)
  setsizes.temp = matrix(0, ncol=3, nrow=5) 
  kk =  5
  f.ase = list()
  for(j in 1:5){
    inds = pData$ID[c(((1:4)+4*(j-1)), ((1:4)+4*(j-1)) + 20, ((1:4)+4*(j-1)) + 40)]
    ki = ((1:4)+4*(j-1))
    nj=length(cor.t2t3[[j]])
    f.ase[[j]] = (rowSums(exprs.all.filt.max[[j]][(1:nj)*2 -1 ,7:18] >= kk) == 12 ) | (rowSums(exprs.all.filt.max[[j]][(1:nj)*2  ,7:18] >= kk) == 12 )
    
    m = match( ratios.max.genes[[j]][f.ase[[j]],4], scaffoldsX)
    f.nx = !is.na(m) 
    #f.ase[[j]] = f.nx  
    # &  f.ase[[j]] # no X 
    q1 = which(f.nx) 
    q2 = q1 
    q3 = q1
    #q1 = which( (rowSums(ratios.max[[j]][f.ase[[j]],-(1:4)] > k | ratios.max[[j]][f.ase[[j]],-(1:4)] < 1-k ) >= 1) )
    #q2 = which((rowSums(ratios.max[[j]][f.ase[[j]],-(5:8)] > k | ratios.max[[j]][f.ase[[j]],-(5:8)] < 1-k ) >= 1)  )
    #q3 = which((rowSums(ratios.max[[j]][f.ase[[j]],-(9:12)] > k | ratios.max[[j]][f.ase[[j]],-(9:12)] < 1-k ) >= 1) )
    
    t1 = prediction_gene_scores( X.tr1[[j]][f.ase[[j]],], X.tt1[[j]][f.ase[[j]],], q1) 
    t2 = prediction_gene_scores( X.tr2[[j]][f.ase[[j]],], X.tt2[[j]][f.ase[[j]],], q2) 
    t3 = prediction_gene_scores( X.tr3[[j]][f.ase[[j]],], X.tt3[[j]][f.ase[[j]],], q3) 
    
    pred.temp[j,] = c(t1[[2]],t2[[2]],t3[[2]])
    pred.ind.temp[j,match(c(t1[[3]],t2[[3]],t3[[3]]), inds)] = 1
    pval.temp[j,] = pval.raw[match(pred.temp[j,], pval.raw[,1]),2]
    setsizes.temp[j,] = c(length(q1), length(q2), length(q3)) 
  }

save(pred.temp, pred.ind.temp, pval.temp,  setsizes.temp, file="ase_identity_xscaff.Rdata")
 return(mean(pred.temp[2:4,]) )
} 
```



## Feature genes  
```{r}
nr = 1:11
nr2 = (10:20)/20
xlab="ratio threshold"

pred.ase.list = list() 
pval.ase.list = list()
pred.indase.list = list() 
setsizes.ase.list = list()
i = 1 
for(k in nr2 ){ 
  #f.ase = lapply(1:5, function(j) (which(rowSums(exprs.all.filt[[j]][,7:18] >= k) == 12 )) ) 
  pred.temp = matrix(0, ncol=3, nrow=5)
  pval.temp = matrix(0, ncol=3, nrow=5)
  pred.ind.temp = matrix(0, ncol=3*4, nrow=5)
  setsizes.temp = matrix(0, ncol=3, nrow=5) 
  kk =  5
  f.ase = list()
  for(j in 1:5){
    inds = pData$ID[c(((1:4)+4*(j-1)), ((1:4)+4*(j-1)) + 20, ((1:4)+4*(j-1)) + 40)]
    ki = ((1:4)+4*(j-1))
    nj=length(cor.t2t3[[j]])
    f.ase[[j]] = (rowSums(exprs.all.filt.max[[j]][(1:nj)*2 -1 ,7:18] >= kk) == 12 ) | (rowSums(exprs.all.filt.max[[j]][(1:nj)*2  ,7:18] >= kk) == 12 )
    
    m = match( ratios.max.genes[[j]][,4], scaffoldsX )
    f.nx = is.na(m) 
    f.ase[[j]] = f.nx  &  f.ase[[j]] # no X 
    
    q1 = which( (rowSums(ratios.max[[j]][f.ase[[j]],-(1:4)] > k | ratios.max[[j]][f.ase[[j]],-(1:4)] < 1-k ) >= 1) )
    q2 = which((rowSums(ratios.max[[j]][f.ase[[j]],-(5:8)] > k | ratios.max[[j]][f.ase[[j]],-(5:8)] < 1-k ) >= 1)  )
    q3 = which((rowSums(ratios.max[[j]][f.ase[[j]],-(9:12)] > k | ratios.max[[j]][f.ase[[j]],-(9:12)] < 1-k ) >= 1) )
    
    t1 = prediction_gene_scores( X.tr1[[j]][f.ase[[j]],], X.tt1[[j]][f.ase[[j]],], q1) 
    t2 = prediction_gene_scores( X.tr2[[j]][f.ase[[j]],], X.tt2[[j]][f.ase[[j]],], q2) 
    t3 = prediction_gene_scores( X.tr3[[j]][f.ase[[j]],], X.tt3[[j]][f.ase[[j]],], q3) 
    
    pred.temp[j,] = c(t1[[2]],t2[[2]],t3[[2]])
    pred.ind.temp[j,match(c(t1[[3]],t2[[3]],t3[[3]]), inds)] = 1
    pval.temp[j,] = pval.raw[match(pred.temp[j,], pval.raw[,1]),2]
    setsizes.temp[j,] = c(length(q1), length(q2), length(q3)) 
  }
  pred.ase.list[[i]] = pred.temp
  pval.ase.list[[i]] = pval.temp
  pred.indase.list[[i]] = pred.ind.temp
  setsizes.ase.list[[i]] = setsizes.temp
  
  
  i = i + 1
}



cd = sapply(nr, function(i) mean(pred.ase.list[[i]]))
bc = sapply(nr, function(i) colMeans(pred.ase.list[[i]]))
ab = sapply(nr, function(i) rowMeans(pred.ase.list[[i]]))
de = sapply(nr, function(i) (pred.ase.list[[i]]))
de.pval = sapply(nr, function(i) (pval.ase.list[[i]]))

bc.pval = matrix(pval.csum[match(bc, pval.csum[,1]),2], nrow=3)
ab.pval = matrix(pval.rsum[match(ab, pval.rsum[,1]),2], nrow=5)
cd.pval = pval.sum[match(cd, pval.sum[,1]),2]

plot( (nr2), cd, pch=19, ylim=c(0,4), main="Gene set aggregate score averaged", ylab="Score", xlab=xlab)
nr22 = sapply(nr, function(i) mean(setsizes.ase.list[[i]]))

cbind(nr, nr2, nr22, cd,t(ab))

save(cd , bc , ab , de , de.pval, bc.pval, ab.pval, cd.pval, nr, nr2, nr22,  pred.ase.list, pval.ase.list, pred.indase.list,  setsizes.ase.list, file="ase_identity_noX.Rdata")   

```

