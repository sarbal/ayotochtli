args <- commandArgs(trailingOnly = TRUE)
n = as.numeric(args[1])

load("cpm_test_train.Rdata")

fracs = 0:100/100 
filt = f.zz
rand.dist = list() 
fi = 1
for(frac in fracs[22:51]){
  
  rand.scores.list= list()
  print(frac)
  for (r in 1:100){    
    X.temp = lapply(1:15, function(i) shuffle_rows(X.cpm.all[filt,(1:4) + 4*(i-1) ]))
    X.sim = do.call(cbind, X.temp)
    
    sampled = sample(sum(filt),n) 
    tempmat = matrix(1:4, ncol= 60, nrow=n, by=T)
    
    X.sim[sampled,] = X.sim[sampled,] + (X.sim[sampled,] * tempmat * frac )
    
    X.rank = lapply(1:15, function(i)   t(apply(X.sim[,(1:4) + 4*(i-1) ],  1, rank )))   
    X.rank2 = do.call(cbind, X.rank)
    
    X.tt1 = X.rank2[,(1:20)]
    X.tt2 = X.rank2[,(21:40)]
    X.tt3 = X.rank2[,(41:60)]
    
    X.t1 = X.sim[,-(1:20)]
    X.t2 = X.sim[,-(21:40)]
    X.t3 = X.sim[,-(41:60)]
    
    X.t1.q = (X.t1[,1:20] + X.t1[,21:40] )/2
    X.t2.q = (X.t2[,1:20] + X.t2[,21:40] )/2
    X.t3.q = (X.t3[,1:20] + X.t3[,21:40] )/2
    
    X.t1.r = lapply(1:5, function(i) t(apply(X.t1.q[,(1:4) + 4*(i-1) ],  1, rank) ))
    X.t2.r = lapply(1:5, function(i) t(apply(X.t2.q[,(1:4) + 4*(i-1) ],  1, rank) ))
    X.t3.r = lapply(1:5, function(i) t(apply(X.t3.q[,(1:4) + 4*(i-1) ],  1, rank) ))
    
    X.tr1 = do.call(cbind,X.t1.r )
    X.tr2 = do.call(cbind,X.t2.r )
    X.tr3 = do.call(cbind,X.t3.r )
    
    
    pred.temp = matrix(0, ncol=3, nrow=5)
    for(j in 1:5){
      inds = pData$ID[c(((1:4)+4*(j-1)), ((1:4)+4*(j-1)) + 20, ((1:4)+4*(j-1)) + 40)]
      ki = ((1:4)+4*(j-1))
      q1 = which(rowSums(X.tt3[,ki] == X.tt2[,ki]) == 4 )
      q2 = which(rowSums(X.tt1[,ki] == X.tt3[,ki]) == 4 )
      q3 = which(rowSums(X.tt1[,ki] == X.tt2[,ki]) == 4 )
      
      t1 = prediction_gene_scores( X.tr1[,ki], X.tt1[,ki], q1) 
      t2 = prediction_gene_scores( X.tr2[,ki], X.tt2[,ki], q2) 
      t3 = prediction_gene_scores( X.tr3[,ki], X.tt3[,ki], q3) 
      pred.temp[j,] = c(t1[[2]],t2[[2]],t3[[2]])
    }
    #mean(pred.temp)
    rand.scores.list[[r]] = pred.temp 
  } 
  rand.dist[[fi]] = rand.scores.list
  fi = fi + 1 
}

save(rand.dist,  fracs, n, file=paste0("rand.distrb.",n,".Rdata" ))


