prediction_gene_scores <- function(training, testing, geneset){
  t1 = 1:dim(training)[2]
  t2 = 1:dim(testing)[2]
  if( length(geneset ) < 2){
    temp12 = matrix(0, ncol=length(t1), nrow=length(t1))
    
    q.sums.cross =  list(temp12, 0, "")
    return(q.sums.cross)
  }
  temp12 = sapply(t1, function(i) colSums(training[geneset,i] == testing[geneset,t2] )  )
  id_score = count_identity_mod(temp12)
  q.sums.cross =  list(temp12, id_score[[1]], id_score[[2]]  )
  return(q.sums.cross)
}


count_identity_mod <- function( temp1) {
  score = 0
  id.list=list()
  colnames(temp1) = rownames(temp1)
  temp1[is.na(temp1)] = 0
  if(sum(temp1) == 0) {  return (list(0, "")) }
  
  shuff = sample( dim(temp1)[1])  
  temp1 = temp1[shuff,shuff]
  for(i in 1:(dim(temp1)[1]-1)){
    pred1 = t(which(temp1== temp1[which.max(temp1)], arr.ind=T))
    p1 = rownames(temp1)[pred1[1]]
    p2 = colnames(temp1)[pred1[2]]
    
    if( p1 == p2){
      score = score +1
      id.list[[i]] = p1
    }
    p1left = rownames(temp1)[-pred1[1]]
    p2left = colnames(temp1)[-pred1[2]]
    temp1 = temp1[-pred1[1], -pred1[2]]
    
  }
  if(p1left == p2left ) {
    score = score +1
    id.list[[i+1]] = p1left
  }
  id.list = unlist(id.list)
  return (list(score, id.list))
}

convolve_x <- function(x,xc,n){
  xcc =  convolve(xc, rev(x), type="open")
  if (n ==1){
    return(xcc)
  } else {
    return(convolve_x(x, xcc, n-1) ) 
  }
}


count_identities_analytic <- function(n){
  A = permutations(n=n,r=n, v=1:n) 
  scores = matrix(0, ncol=1,nrow=dim(A)[1])
  for( j in 1:dim(A)[1] ){ 
    pred1 = 1:n
    pred2 = A[j,]
    score = 0
    temp1 = diag(n) * 0 
    colnames(temp1) = 1:n
    rownames(temp1) = 1:n
    
    #    pred1 = rbind(pred1,pred2)
    for(i in 1:n){ 
      temp1[pred1[i],pred2[i] ] = 1 
    }
    
    
    for(i in 1:(n-1)) {
      pred1 = t(which(temp1== temp1[which.max(temp1)], arr.ind=T))
      
      # predictor
      p1 = rownames(temp1)[pred1[1]]
      # predictions 
      p2 = colnames(temp1)[pred1[2]]
      
      if( p1 == p2){
        score = score +1
      }
      p1left = rownames(temp1)[-pred1[1]]
      p2left = colnames(temp1)[-pred1[2]]
      temp1 = temp1[-pred1[1], -pred1[2]]
      
      
    }
    if(p1left == p2left ) { 
      score = score +1
    }  
    scores[j,1] = score
  }
  return(scores)
}





prediction_gene_scores_mod <- function(training, testing, geneset){
  t1 = 1:dim(training)[2]
  t2 = 1:dim(testing)[2]
  if ( length(geneset) < 1 ){ 
    temp12 = matrix(0, ncol=length(t1), nrow=length(t1))
    q.sums.cross =  list(temp12, 0, "")
    return(q.sums.cross)
  }
  if( length(geneset ) == 1){
   temp12 = sapply(t1, function(i) 1*(training[geneset,i] == testing[geneset,t2] )  )
    
  } else { 
  temp12 = sapply(t1, function(i) colSums(training[geneset,i] == testing[geneset,t2] )  )
 } 
  id_score = count_identity_mod(temp12)
  q.sums.cross =  list(temp12, id_score[[1]], id_score[[2]]  )
  return(q.sums.cross)
  
}

prediction_gene_scores_indiv <- function(training, testing, geneset){
  t1 = 1:dim(training)[2]
  t2 = 1:dim(testing)[2]
 temp = list() 
  for(j in geneset) {
   temp12 = sapply(t1, function(i) 1*(training[j,i] == testing[j,t2] )  )
   id_score = count_identity_mod(temp12)
   temp[[j]] = id_score[[1]]
  } 
   return(unlist(temp) ) 
 
  
}

