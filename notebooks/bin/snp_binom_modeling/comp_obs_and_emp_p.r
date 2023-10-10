require(ggplot2)
require(ggsci)
require(viridis)
require(tidyverse)
require(dplyr)
rdf <- readRDS('simple_independent_result.rds')
sdf <- readRDS('simple_t0_result.rds')
# sig_df <- readRDS('simple_independent_sig_summary.rds')
# sig_t0 <- readRDS('simple_t0_sig_summary.rds')
all_df <- readRDS('simple_independent_all_summary.rds')
all_t0 <- readRDS('simple_t0_all_summary.rds')
blacklist <- readRDS("black_list_homo_snps.rds")
load("data/Xscaffolds.Rdata")

add.two.xscaffold <- function() {
    scaffoldsX.prop <- c(scaffoldsX.prop, "JH573670", "JH583104")
    scaffoldsX.prop2 <- c(scaffoldsX.prop2, "JH573670", "JH583104")
    scaffoldsX.prop3 <- c(scaffoldsX.prop3, "JH573670", "JH583104")
}
add.two.xscaffold()
fix.prob.one.sample <- function(data) {
    data[data[,'pvalue'] == 1,'pvalue'] = 0.5
    data[data[,'pvalue_fixed'] == 1, 'pvalue_fixed'] = 0.5
    return(data)
}
fix.prob.one <- function(df) {
    for (ident in 1:length(df)) {
        for (time in 1:length(df[[ident]])) {
            tdata <- as.data.frame(df[[ident]][[time]])
            print(head(tdata))
            print(colnames(tdata))
            df[[ident]][[time]] <- fix.prob.one.sample(tdata)
        }
    }
    return(df)
}
# rdf <- fix.prob.one(rdf)
# sdf <- fix.prob.one(sdf)

print(head(rdf[[1]][[1]]))
print(head(sdf[[1]][[1]]))
# print(head(sig_df))
# print(head(sig_t0))
print(head(all_df))

read.emp.result <- function(i) {
    df <- readRDS(paste0('simple_independent_result_emp_', i, '.rds'))
    df <- fix.prob.one.sample(df)
    return(df)
}
remove.homozygous.snps <- function(df, quad_bl, fixed, verbose=TRUE) {
    min_c = 10
    if (verbose) {
        print('remove by homozygous')
        print(dim(df))
    }
    df <- df[unlist(sapply(df[,'gene'], function(x){return(!any(x == quad_bl))})), ]
    if (verbose) print(dim(df))
    # p > 0.9 and min cov for minor allele 
    if (fixed) {
        flag = (df[,'p_fixed'] > 0.9 | df[,'p_fixed'] < 0.1)
    } else {
        flag = (df[,'p'] > 0.9 | df[,'p'] < 0.1)
    }
    remove <- c()
    for (g in unique(df[flag,'gene'])) {
        index = which(df[,'gene'] == g)
        index = index[df[index,'cov'] > 0]
        actual_count = df[index, 'x']*df[index, 'cov']
        # print(actual_count)
        if (any(!is.finite(actual_count))) {
            print(df[index,])
        }
        if (max(actual_count) < min_c || max(df[index, 'cov']-actual_count) < min_c)
            remove <- c(remove, g)
    }
    #print(c('blacklist', remove))
    df <- df[!(df[,'gene'] %in% remove),]
    if (verbose) print(dim(df))
    return(df)
}

plot.obs.emp.stat.remove.x <- function(df, header, ident_index, fixed=FALSE) {       
    require(dplyr)
    print(header)
    ident = floor((ident_index-1)/3)+1
    # print(c(ident, ident_index))
    quad_bl = blacklist[[ident]]
    temp = all_df[(all_df[,'ident'] == ident) & (all_df[,'time'] == 1),]
    print(head(all_df))
    print(head(temp))
    stopifnot(dim(temp)[1] > 0)
    # temp = rdf[[ident]][[(ident_index-1)%%3+1]]
    remained <- temp[!(temp[,'chrm'] %in% scaffoldsX.prop),'gene']
    print(scaffoldsX.prop)
    print(length(remained))
    print(head(remained))
    print(head(df))
    tdf <- df[df[,'gene'] %in% remained,]
    print('removed by x scaffolds')
    print(ident_index)
    print(dim(df))
    print(dim(tdf))
    plot.obs.emp.stat(tdf, paste0('plot_', header, '_remx'), paste0(header, '_remx', '.tsv'), quad_bl, fixed=fixed)
}

summarize.df <- function(df) {
    tdf <- df %>% group_by(gene) %>% summarise(
        mean_cov = mean(cov),
        # cor      = min(cor),
        # scor     = min(scor),
        prob     = -sum(prob, na.rm=TRUE),
        pvalue   = -sum(pvalue+log(2), na.rm=TRUE),
        prob_fixed   = -sum(prob_fixed, na.rm=TRUE),
        pvalue_fixed = -sum(pvalue_fixed+log(2), na.rm=TRUE),
        p = min(p),
        p_fixed = min(p_fixed),
        fisher_obs = min(fisher_obs),
        fisher_emp = min(fisher_emp),
        fisher_obs_fixed = min(fisher_obs_fixed),
        fisher_emp_fixed = min(fisher_emp_fixed),
        count = n()
    )
    tdf <- as.data.frame(tdf)
    return(tdf)
}
compute.fisher.fdr <- function(tdf, header, x_label, y_label, plot=TRUE) {
    emp_obs <- do.call(rbind, lapply(1:dim(tdf)[1], function(x) {
        return(c(-log10(max(0.0001, 1.-tdf[x,x_label]/10000)), min(-log(0.0001), -(pchisq(tdf[x,y_label], 2*(tdf[x,'count'])-1, lower.tail = FALSE, log.p = TRUE)))/log(10)))
    }))
    emp_obs <- data.frame(emp_obs)
    colnames(emp_obs) <- c('emp_p', 'chi_p')
    if (plot) {
        g <- ggplot(emp_obs, aes(x=emp_p, y=chi_p))+geom_point()+geom_density2d()+
        stat_density2d(aes(alpha=..level.., fill=..level..), size=2, 
            bins=10, geom="polygon") + scale_fill_gradient(low = "yellow", high = "red") + theme_bw()
        png(paste0(header, '.png'))
        plot(g)
        dev.off()
        g <- ggplot(emp_obs, aes(x=emp_p, y=chi_p))+geom_point()+geom_density2d()+
        stat_density2d(aes(alpha=..level.., fill=..level..), size=2, 
            bins=100, geom="polygon") + scale_fill_gradient(low = "yellow", high = "red") + theme_bw()+xlim(0, 0.1)
        png(paste0(header, '_around_zero.png'))
        plot(g)
        dev.off()
    }
    emp_obs <- cbind(emp_obs, emp_p_adj=p.adjust(10**(-emp_obs[,1]), method='fdr'))
    emp_obs <- cbind(emp_obs, chi_p_adj=p.adjust(10**(-emp_obs[,2]), method='fdr'))
    if (plot) {
        g <- ggplot(emp_obs, aes(x=emp_p_adj, y=chi_p_adj))+geom_point()+geom_density2d()+scale_x_continuous(trans='log10') +
    scale_y_continuous(trans='log10') + stat_density2d(aes(alpha=..level.., fill=..level..), size=2, 
            bins=10, geom="polygon") + scale_fill_gradient(low = "yellow", high = "red") + theme_bw()
        png(paste0(header, '_adj.png'))
        plot(g)
        dev.off()
    }
    return(emp_obs)
}

plot.obs.emp.stat <- function(df, header, output, quad_bl=c(), fixed=FALSE) {       
    require(ggplot2)
    require(viridis)
    require(dplyr)
    require(reshape2)
    x_label <- 'fisher_emp'
    y_label <- 'fisher_obs'
    if (fixed) {
        x_label <- paste0(x_label, '_fixed')
        y_label <- paste0(y_label, '_fixed')
    }
    df <- remove.homozygous.snps(df, quad_bl, fixed)
    tdf <- summarize.df(df)
    emp_obs <- compute.fisher.fdr(tdf, header, x_label, y_label,  plot=FALSE)
    write.table(cbind(tdf, emp_obs), file=output, sep='\t', quote=F)
    mdf <- melt(emp_obs[,c('emp_p', 'chi_p')])
    print(head(mdf))
    print(dim(subset(mdf, mdf[,'variable'] == 'emp_p')))
    print(dim(subset(mdf, mdf[,'variable'] == 'chi_p')))
    mdf[,'value'] = 10**(-mdf[,'value'])
    g <- ggplot(mdf, aes(x=value, fill=variable))+geom_histogram(alpha=0.5, position='identity') + theme_bw()
    png(paste0(header, '_hist.png'))
    plot(g)
    dev.off()
    g <- ggplot(mdf, aes(x=value, fill=variable))+geom_histogram(alpha=0.5, position='identity', bins=50) + theme_bw()+xlim(0, 0.005)
    png(paste0(header, '_hist_around_zero.png'))
    plot(g)
    dev.off()
    print(dim(tdf))
    print(dim(emp_obs))
}

for (i in 1:15) {
    df <- read.emp.result(i)
    plot.obs.emp.stat(df, paste0('plot_fisher_stat_obs_emp_', i), paste0('fisher_stat_obs_emp_', i, '.tsv'), blacklist[[(i-1)/3+1]])
    plot.obs.emp.stat(df, paste0('plot_fisher_stat_obs_emp_fixed_', i), paste0('fisher_stat_obs_emp_fixed_', i, '.tsv'), blacklist[[(i-1)/3+1]], fixed=TRUE)
    plot.obs.emp.stat.remove.x(df, paste0('fisher_stat_obs_emp_', i), i, fixed=FALSE)
    plot.obs.emp.stat.remove.x(df, paste0('fisher_stat_obs_emp_fixed_', i), i, fixed=TRUE)
}
