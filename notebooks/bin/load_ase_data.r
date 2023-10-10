
load("ase_ratios.test_train.Rdata")
load("exprs_all.Rdata")
load("metadata.Rdata")
load("gene_annotations_v0.95.Rdata")
load("armadillo.helper.Rdata")
# Variables 
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



