require(ggplot2)
require(ggridges)
require(viridis)
require(psych)
require(ggpubr)
require(stringr)
require(viridis)
require(ggsci)
require(dplyr)
require(tidyverse)
require(ComplexHeatmap)
require(reshape2)
require(parallel)

homo.prob <- function(p) {
    return(p*p)
}
hetero.prob <- function(p) {
    return(p*(1-p)*2)
}

# use the average of left and right edges of each window
MEAN = TRUE
# maximum cell number
max_cell <- 1000000
# p of binomial distribution
p <- 0.5
# params for beta distribution
params_a <- c(1, 2, 4, 10, 15, 30, 1, 1, 1, 4, 10, 30)
params_b <- c(1, 2, 4, 10, 15, 30, 4, 10, 30, 1, 1, 1)
step <- 0.000001
start <- 0
all_df <- NULL

for (count in 1:length(params_a)) {
    if (count < 10) next # compute the models 11, 12, 13
    temp <- NULL
    mcount <- 0
    mutation_rate = seq(1, 5, 0.1
    for (i in 1:length(mutation_rate)) {
        mut_p = mutation_rate[i]
        sig_ratio <- 10**(-mut_p)
        ccount <- 0
        for (cnum in c(2, 4, 8, 16, 32, 64, 128, 256, max_cell)) {
	    # compute integral of binomial (diff of pbinom) and cumulative beta (average or one point)
            if (MEAN) {
                x = seq(start, 1.0, step)
                y = unlist(mclapply(x, function(q){
                    right = pbeta(q, params_a[count], params_b[count], ncp = 0, lower.tail = TRUE, log.p = FALSE)
                    if (q == 0.) {
                        left = right
                        left_b = 0
                    } else {
                        left = pbeta(q-step, params_a[count], params_b[count], ncp = 0, lower.tail = TRUE, log.p = FALSE)
                        left_b = pbinom((q-step)*cnum, size=cnum, prob=p, lower.tail=TRUE)
                    }
                    stopifnot(right >= left)
                    right_b = pbinom(q*cnum, size=cnum, prob=p, lower.tail=TRUE)
                    if (q == 0) left_b = 0
                    stopifnot(right_b >= left_b)
                    return((right+left)/2*(right_b-left_b))
                }, mc.cores=20))
            } else {
                x = seq(start, 1.0, step)
                y = unlist(mclapply(x, function(q){
                    right = pbeta(q, params_a[count], params_b[count], ncp = 0, lower.tail = TRUE, log.p = FALSE)
                    right_b = pbinom(q*cnum, size=cnum, prob=p, lower.tail=TRUE)
                    if (q == 0) left_b = 0
                    else left_b = pbinom((q-step)*cnum, size=cnum, prob=p, lower.tail=TRUE)
                    return(right*(right_b-left_b))
                }, mc.cores=20))
            }
            stopifnot(sum(y) >= 0)
            one_data = data.frame(a=params_a[count], b=params_b[count], total=sum(y)*hetero.prob(sig_ratio)+homo.prob(sig_ratio), 
                                    both=sum(y)*hetero.prob(sig_ratio)/(homo.prob(sig_ratio)+sum(y)*hetero.prob(sig_ratio)), 
                                    hetero=sum(y)*hetero.prob(sig_ratio), count=count, cnum=paste0(ccount, '_', cnum), mutation_rate=paste0(mcount, '_', round(10**(-mut_p), 4)))
            if (is.null(temp)) {
                temp = one_data
            } else {
                temp = rbind(temp, one_data)
            }
            ccount = ccount+1
        }
        mcount = mcount+1
    }
    for (prob in c('both', 'hetero', 'total')) {
        ttemp = temp[,c('mutation_rate', 'cnum', prob)]
        tail = ''
        if (MEAN) tail = '_mean'
        colnames(ttemp)[3] = 'prob'
        print(head(temp))
        mat = dcast(ttemp, mutation_rate ~ cnum)
        print(head(mat))
        rownames(mat) = mat[,'mutation_rate']
        mat = mat[,2:dim(mat)[2]]
        mat = as.matrix(mat)
        heatmap = Heatmap(mat, col=viridis(100), cluster_rows=FALSE, cluster_columns=FALSE)
        pdf(paste0('plot_pdist_', prob, '_', count, tail, '.pdf'))
        draw(heatmap)
        dev.off()
        mat[mat == 0] = 1e-20
        mat = log10(mat)
        heatmap = Heatmap(mat, col=viridis(100), cluster_rows=FALSE, cluster_columns=FALSE)
        pdf(paste0('plot_pdist_', prob, '_', count, tail, '_log.pdf'))
        draw(heatmap)
        dev.off()
    }
}
