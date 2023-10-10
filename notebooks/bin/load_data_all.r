load("cpm_strand_comb.Rdata")
load("counts_strand_comb.Rdata")
load("gene_annotations_v0_95_mod.Rdata")
load("armadillo_helper.Rdata")
load("armadillo_hem.Rdata")
load("cpm_test_train.Rdata")
dataset = "cpm"




# Variables 
N = dim(X.cpm.all)[1]
n_quads = 5
n_times = 3 
n_qt = n_quads * n_times
n_q = 4 
n_samp = n_quads * n_q
r_samp = 1:n_samp 

# Labels 
quad = as.numeric(pData$Quad)
sex = as.numeric(pData$Sex)
lane = as.numeric(pData$Batches)
tlab = c("t1", "t2", "t3")
quads = unique(substr(pData$ID, 1,4))

labels = paste(sapply(1:n_quads, function(i) rep(quads[i],n_times)), tlab  ) 
sexlabels = unlist(lapply(1:n_quads, function(i) rep(unique(cbind(pData$Quad, pData$Sex) )[i,2],n_times)) ) 
timelabels = rep(1:n_times, n_quads)
quadlabels = unlist(lapply(1:n_quads, function(i) rep(unique(cbind(pData$Quad, pData$Sex) )[i,1],n_times))   ) 


# Z-scores per quad per timepoint 
X.zscores = lapply(1:n_qt, function(i)  t(scale( t(X.cpm.all[, (1:4) + 4*(i-1)  ]))))   
X.zscores = do.call(cbind, X.zscores)


# Ranks per quad per timepoint (across quads)
X.rank = lapply(1:n_qt, function(i)   t(apply(X.cpm.all[,(1:4) + 4*(i-1) ],  1, rank )))   
X.rank2 = do.call(cbind, X.rank)

# Gene filters 
f.z = sapply(1:n_qt, function(i) rowSums(X.cpm.all[,1:4+ 4*(i-1)]) > 0  )
f.zz = rowSums(f.z) == n_qt 
f.zzz = (rowSums(X.cpm.all == 0 )  == 0  )
