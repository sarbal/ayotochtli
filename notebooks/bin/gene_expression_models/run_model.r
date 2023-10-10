args <- commandArgs(trailingOnly = TRUE)
    j = as.numeric(args[1])

load("run_model.Rdata")
k = 0.8 
sds.data = sample( all.sds , Nmax)
rand.ngenes = list() 
ncells = c(1,2,4,8,16,32,64,100,200,300,500,1000)
fracs = 1:50/100
fracs = fracs/10

nsamp = round(Nmax/2 * fracs[j])
rand.N = list()
for (N in ncells){ 
  print(N)
  rand.mod = list() 
  for( r in 1:1000){ 
        ratio.sim = matrix(rnorm(8*Nmax,mean=0.5,sd=sqrt(0.25/N)), ncol=8, nrow=Nmax)
        ratio.sim[ratio.sim>1] = 1
        ratio.sim[ratio.sim<0] = 0
        
        #hist(ratio.sim, xlim=c(0,1))
        
        s1 = sample(Nmax,nsamp)
        s2 = sample(Nmax,nsamp)
        
        X1 = ratio.sim[s1,1:4] 
        X2 = ratio.sim[s2,1:4]  
        
        #X1 = ratio.sim[s1,1:4] * c(1.1,1.2,1.3, 1.4)
        #X2 = ratio.sim[s2,1:4] * c(0.9,0.8,0.7, 0.6)
        
        # Xi1 = sapply(1:4, function(i) apply( as.matrix(X1[,i]), 1, pick_dist ,  sqrt(0.25/N)))
        # Xi2 = sapply(1:4, function(i) apply( as.matrix(X2[,i]), 1, pick_dist ,  sqrt(0.25/N)))
        
        Xi1 = sapply(1:4, function(i) apply( as.matrix(X1[,i]), 1, pick_dist ,  sds.data[s1]))
        Xi2 = sapply(1:4, function(i) apply( as.matrix(X2[,i]), 1, pick_dist ,  sds.data[s2]))

        Xii = rbind(cbind(X1,Xi1), cbind(X2,Xi2) )
        
        ratio.sim[c(s1,s2),] = Xii
        # hist(ratio.sim, xlim=c(0,1))
        
        X.tr1 = t(apply(ratio.sim[,1:4], 1, rank))
        X.tt1 = t(apply(ratio.sim[,5:8], 1, rank))
        
        colnames(X.tr1) = 1:4
        colnames(X.tt1) = 1:4
        q1 = which(rowSums(ratio.sim >= k | ratio.sim <= 1-k ) > 0)
        t1 = prediction_gene_scores( X.tr1, X.tt1,q1)
        rand.mod[[r]] = t1[[2]]
  }

  rand.N[[N]] = rand.mod
} 



save(rand.N, ncells, fracs, j, file=paste0("rand.fracs2.",j,".Rdata" ))


