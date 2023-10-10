require(ggplot2)
require(ggsci)
require(viridis)
require(tidyverse)
require(parallel)
require(dplyr)
rdf <- readRDS('simple_independent_result.rds')
sdf <- readRDS('simple_t0_result.rds')
sig_df <- readRDS('simple_independent_sig_summary.rds')
sig_t0 <- readRDS('simple_t0_sig_summary.rds')
all_df <- readRDS('simple_independent_all_summary.rds')
all_t0 <- readRDS('simple_t0_all_summary.rds')
load("data/Xscaffolds.Rdata")
add.two.xscaffold <- function() {
    scaffoldsX.prop <- c(scaffoldsX.prop, "JH573670", "JH583104")
    scaffoldsX.prop2 <- c(scaffoldsX.prop2, "JH573670", "JH583104")
    scaffoldsX.prop3 <- c(scaffoldsX.prop3, "JH573670", "JH583104")
}
add.two.xscaffold()
print(head(rdf[[1]][[1]]))
print(head(sdf[[1]][[1]]))
print(head(sig_df))
print(head(sig_t0))
print(head(all_df))

compute.emp.pvalue <- function(source, header, ind) {
    cores = 50
    for (target in c('', '_fixed')) {
        emp.pvalue <- do.call(rbind, mclapply(unique(source[,'gene']), function(x) {
            if (x %%10000 == 1) print(x)
            temp <- source[source[,'gene'] == x,]
            total_cov = sum(temp[,'cov'])
            if (target == '')
                obs_p = temp[1,'p']
            else
                obs_p = temp[1,'p_fixed']
            obs_stat = -2.*sum(unlist(sapply(1:dim(temp)[1], function(j) {
                return(log(2)+pbinom(temp[j,'x']*temp[j,'cov'], temp[j, 'cov'], obs_p, log.p=TRUE, lower.tail=(temp[j,'x']*temp[j,'cov'] < obs_p)))
            })))
            samp_stat = sapply(1:10000, function(i) {
                sample_p = unlist(sapply(1:dim(temp)[1], function(j) {
                    return(rbinom(1, temp[j,'cov'], obs_p))
                }))
                if (target == '')
                    emp_p = sum(sample_p)/total_cov
                else
                    emp_p = mean(sample_p/temp[,'cov'], na.rm=TRUE)
                emp_stat = -2.*sum(unlist(sapply(1:dim(temp)[1], function(j) {
                    return(log(2)+pbinom(sample_p[j], temp[j,'cov'], emp_p, log.p=TRUE, lower.tail=(sample_p[j] < emp_p)))
                })))
                return(emp_stat)
            })
            if (x%%10000 == 1) {
                png(paste0('comp_theory_obs_emp_', header, '_', x, target, '.png'))
                chi_sample = rchisq(1000, dim(temp)[1]*2-1)
                total_sample = c(chi_sample[is.finite(chi_sample)], samp_stat[is.finite(samp_stat)])
                hist(samp_stat, col=rgb(1, 0, 0, 1/4), freq=FALSE, xlim=c(min(total_sample), max(total_sample)))
                hist(chi_sample, col=rgb(0, 0, 1, 1/4), freq=FALSE, add=TRUE)
                abline(v=obs_stat, col='red')
                dev.off()
            }
            rank = length(which(samp_stat <= obs_stat))
            return(cbind(fisher_obs=rep(obs_stat, dim(temp)[1]), fisher_emp=rep(rank, dim(temp)[1])))
        }, mc.cores=cores))
        print(dim(emp.pvalue))
        print(dim(source))
        colnames(emp.pvalue) = c(paste0('fisher_obs', target), paste0('fisher_emp', target))
        source <- cbind(source, emp.pvalue)
    }
    saveRDS(source, file=paste0('simple_independent_result_emp_', ind, '.rds'))
    return(source)
}

compute.emp.sim <- function(all_time_points, header, ind) {
    require(lava)
    cores = 50
    num = 100
    for (target in c('', '_fixed')) {
        result = list()
        for (i in 1:3) {
            print(c(ind, i))
            print(dim(all_time_points[[i]]))
            source = all_time_points[[i]]
            # print(length(all_time_points))
            genes = unique(all_time_points[[i]][,'gene'])
            emp.pvalue <- mclapply(genes, function(x) {
                if (x %%10000 == 1) print(x)
                temp <- source[source[,'gene'] == x,]
                total_cov = sum(temp[,'cov'])
                if (target == '')
                    obs_p = temp[1,'p']
                else
                    obs_p = temp[1,'p_fixed']
                obs_stat = -2.*sum(unlist(sapply(1:dim(temp)[1], function(j) {
                    return(log(2)+pbinom(temp[j,'x']*temp[j,'cov'], temp[j, 'cov'], obs_p, log.p=TRUE, lower.tail=(temp[j,'x']*temp[j,'cov'] < obs_p)))
                })))
                samp_stat = do.call(rbind, lapply(1:num, function(i) {
                    sample_p = unlist(sapply(1:dim(temp)[1], function(j) {
                        return(rbinom(1, temp[j,'cov'], obs_p))
                    }))
                    if (target == '')
                        emp_p = sum(sample_p)/total_cov
                    else
                        emp_p = mean(sample_p/temp[,'cov'], na.rm=TRUE)
                    emp_stat = -2.*sum(unlist(sapply(1:dim(temp)[1], function(j) {
                        return(log(2)+pbinom(sample_p[j], temp[j,'cov'], emp_p, log.p=TRUE, lower.tail=(sample_p[j] < emp_p)))
                    })))
                    vec <- c(unlist(sapply(1:dim(temp)[1], function(j){return(sample_p[j]/temp[j,'cov'])})), emp_stat, x, total_cov)
                    if (length(vec) < 7) {
                        return(rep(0, 7))
                    } else return(vec)
                }))
                colnames(samp_stat) = c(1, 2, 3, 4, 'exm_stat', 'gene', 'cov')
                return(samp_stat)
            }, mc.cores=cores)
            for (j in 1:length(emp.pvalue)) {
                rownames(emp.pvalue[[j]]) = paste0(genes[j], '_', 1:dim(emp.pvalue[[j]])[1])
            }
            mat <- do.call(rbind, emp.pvalue)
            result[[i]] <- as.matrix(mat)
        }
        print(length(result))
        print(head(result[[1]]))
        print(head(result[[2]]))
        print("## 0 coverage is excluded")
        for (i in 2:3) {
            print(dim(result[[i]]))
            kept_genes <- rownames(result[[i]])[rownames(result[[i]]) %in% rownames(result[[1]])]
            print(rownames(result[[1]])[!(rownames(result[[1]]) %in% rownames(result[[i]]))])
            print(length(kept_genes))
            result[[i]] = result[[i]][kept_genes,]
        }
        for (i in 2:3) {
            cor_mat = do.call(rbind, mclapply(rownames(result[[i]]), function(x){
                vec = result[[1]][which(rownames(result[[1]]) == x),]
                vec2 = result[[i]][which(rownames(result[[i]]) == x),]
                return(c(cor(vec[1:4], vec2[1:4]), cor(vec[1:4], vec2[1:4], method='spearman')))
            }, mc.cores=cores))
            print(head(cor_mat))
            colnames(cor_mat) <- c('cor', 'scor')
            result[[i]] <- cbind(result[[i]], cor_mat)
        }
        print(head(result[[2]]))
        saveRDS(result, file=paste0('simple_independent_result_emp_sim_', ind, target, '_', num, '.rds'))
        # return(result)
    }
    return(NULL)
}

print(length(rdf))
print(length(rdf[[1]]))
CORE = 1
# new_rdf <- mclapply(1:15, function(x){compute.emp.pvalue(rdf[[as.integer((x-1)/3)+1]][[(x-1)%%3+1]], paste0(as.integer((x-1)/3)+1, '_', (x-1)%%3+1), x)}, mc.cores=CORE)
new_rdf <- mclapply(1:5, function(x){compute.emp.sim(rdf[[as.integer((x-1)/3)+1]], paste0(as.integer((x-1)/3)+1, '_', (x-1)%%3+1), x)}, mc.cores=CORE)
# saveRDS(new_rdf, file='simple_independent_result_emp.rds')