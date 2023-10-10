require(ggplot2)
require(ggsci)
require(viridis)
require(tidyverse)
require(dplyr)
rdf <- readRDS('simple_independent_result.rds')
sdf <- readRDS('simple_t0_result.rds')
sig_df <- readRDS('simple_independent_sig_summary.rds')
sig_t0 <- readRDS('simple_t0_sig_summary.rds')
all_df <- readRDS('simple_independent_all_summary.rds')
all_t0 <- readRDS('simple_t0_all_summary.rds')
load("data/Xscaffolds.Rdata")
blacklist <- readRDS("black_list_homo_snps.rds")
add.two.xscaffold <- function() {
    scaffoldsX.prop <- c(scaffoldsX.prop, "JH573670", "JH583104")
    scaffoldsX.prop2 <- c(scaffoldsX.prop2, "JH573670", "JH583104")
    scaffoldsX.prop3 <- c(scaffoldsX.prop3, "JH573670", "JH583104")
}
add.two.xscaffold()
fix.prob.one <- function(df) {
    for (ident in length(df)) {
        for (time in length(df[[ident]])) {
            tdata <- df[[ident]][[time]]
            tdata[tdata[,'pvalue'] == 1,'pvalue'] = 0.5
            tdata[tdata[,'pvalue_fixed'] == 1, 'pvalue_fixed'] = 0.5
            df[[ident]][[time]] <- tdata
        }
    }
    return(df)
}
# rdf <- fix.prob.one(rdf)
# saveRDS(rdf, 'simple_independent_result_updated.rds')
# sdf <- fix.prob.one(sdf)
# saveRDS(sdf, 'simple_t0_result_updated.rds')

print(head(rdf[[1]][[1]]))
print(head(sdf[[1]][[1]]))
print(head(sig_df))
print(head(sig_t0))
print(head(all_df))

remove.homozygous.snps <- function(df, quad_bl, fixed) {
    min_c = 10
    df <- df[unlist(sapply(df[,'gene'], function(x){return(!any(x == quad_bl))})), ]
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
    return(df)
}

plot.hist.two.correlation <- function(sig_df, all_df, cor, header, label='sig') 
{
    combined <- rbind(all_df, data.frame(cor=sig_df[,cor], snp=rep(label, dim(sig_df)[1])))
    print(dim(combined))
    combined[,2] <- factor(combined[,2])
    print(head(combined))
    png(paste0('hist_', header, '.png'))
    if (cor == 'scor')
        g <- ggplot(combined, aes(x=cor, color=snp, fill=snp))+geom_histogram(aes(y=..density..), alpha=0.5, position='identity', )+scale_color_npg()+theme_bw()
    else
        g <- ggplot(combined, aes(x=cor, color=snp, fill=snp))+geom_histogram(aes(y=..density..), alpha=0.5, position='identity')+scale_color_npg()+theme_bw()
    plot(g)
    dev.off()

}
plot.boxplot.correlation <- function(df, header)
{
    require(ggplot2)
    g <- ggplot(df, aes(x=cor, y=score))+geom_boxplot()+theme_bw()
    # +geom_jitter(shape=16, size=1, position=position_jitter(0.2), alpha=0.3)
    png(paste0('boxplot_', header, '.png'))
    plot(g)
    dev.off()
}
sd.min <- function(x) {
    return(mean(x)+sd(x))
}
sd.max <- function(x) {
    return(mean(x)-sd(x))
}

plot.pn.accum <- function(data, header) {
    data[,'x'] <- rank(data[,'x'])
    data[,'group'] = sapply(data[,'y'], function(x){if(x > 0) return('P'); if (x < 0) return('N'); return(NA)})
    tdf <- as.data.frame(data) %>%
        arrange(x) %>%
        mutate(rolling_sd=rollapply(y, width=width, sd, fill=NA), rolling_mean=rollapply(y, width=width, mean, fill=NA))
    # print(c(min(tdf[,'rolling_mean']-tdf[,'rolling_sd']/2, na.rm=T), max(tdf[,'rolling_mean']+tdf[,'rolling_sd']/2, na.rm=T)))
    g <- ggplot(tdf, aes(x=x, y=y))+geom_point(alpha=0.01)+geom_line(aes(x=x, y=rolling_mean), color='red')+geom_ribbon(aes(ymin=rolling_mean-rolling_sd/2, ymax=rolling_mean+rolling_sd/2), fill="lightgray", color="lightgray", alpha=.8) + theme_bw()
    g <- g + ylim(min(tdf[,'rolling_mean']-tdf[,'rolling_sd']/2, na.rm=T), max(tdf[,'rolling_mean']+tdf[,'rolling_sd']/2, na.rm=T))
    png(paste0('rank_mean_ave_', header, '.png'))
    plot(g)
    dev.off()

}
plot.rank.mean.ave <- function(data, header, rank=FALSE, all=FALSE) {
    require(ggplot2)
    require(zoo)
    width = min(max(2, dim(data)[1]/4/100), 100)
    if (all) width = 500
    data[,'x'] <- rank(data[,'x'])
    if (rank) {
        data[,'y'] <- rank(data[,'y'])
    }
    tdf <- as.data.frame(data) %>%
        arrange(x) %>%
        mutate(rolling_sd=rollapply(y, width=width, sd, fill=NA), rolling_mean=rollapply(y, width=width, mean, fill=NA))
    # print(c(min(tdf[,'rolling_mean']-tdf[,'rolling_sd']/2, na.rm=T), max(tdf[,'rolling_mean']+tdf[,'rolling_sd']/2, na.rm=T)))
    g <- ggplot(tdf, aes(x=x, y=y))+geom_point(alpha=0.01)+geom_line(aes(x=x, y=rolling_mean), color='red')+geom_ribbon(aes(ymin=rolling_mean-rolling_sd/2, ymax=rolling_mean+rolling_sd/2), fill="lightgray", color="lightgray", alpha=.8) + theme_bw()
    g <- g + ylim(min(tdf[,'rolling_mean']-tdf[,'rolling_sd']/2, na.rm=T), max(tdf[,'rolling_mean']+tdf[,'rolling_sd']/2, na.rm=T))
    png(paste0('rank_mean_ave_', header, '.png'))
    plot(g)
    dev.off()
    if (!rank) {
        if (dim(tdf)[1] < 100)
            g <- ggplot(tdf, aes(x=x, y=y))+geom_point()+theme_bw()+scale_fill_continuous(type='viridis')
        else
            g <- ggplot(tdf, aes(x=x, y=y))+geom_hex(bins=50)+theme_bw()+scale_fill_continuous(type='viridis')        
        png(paste0('dens_', header, '.png'))
        plot(g)
        dev.off()
    }
}

plot.p.hist <- function(p_stat, header) {
    require(ggpubr)
    require(ggplot2)
    require(ggsci)
    print(head(p_stat))
    print(unique(p_stat[,'ident']))
    print(unique(p_stat[,'time']))
    p_stat[,'ident'] <- factor(p_stat[,'ident'])
    prev = -1
    for (cov in c(-1, 10, 100, 250, 500, 750, 1000)) {
        if (cov > 0)
            temp <- p_stat[(p_stat[,'cov'] <= cov) & (p_stat[,'cov'] > prev),]
        else
            temp <- p_stat
        prev = cov
        tail = paste0('_', cov)
        if (cov < 0) tail = ''
        for (time in 1:3) {
            print(time)
            ttemp <- temp[temp[,'time'] == time,]
            g <- ggplot(temp, aes(x = p, fill=ident))+geom_histogram(position='dodge', stat='count', alpha=1)+scale_fill_npg()+theme_bw()
            png(paste0(header, '_', time, tail, '.png'))
            plot(g)
            dev.off()
        }
    }
}


plot.basic.statistics <- function(all_snp, source, header) {
    p_stat <- NULL
    for (ident in 1:5) {
        for (time in 1:length(source[[ident]])) {
            print(c(ident, time))
            # print(head(source[[ident]][[time]]))
            df = source[[ident]][[time]] %>% group_by(gene) %>% summarise(
                mean_cov = mean(cov),
                prob     = -sum(prob, na.rm=TRUE),
                pvalue   = -sum(pvalue+log(2), na.rm=TRUE),
                prob_fixed   = -sum(prob_fixed, na.rm=TRUE),
                pvalue_fixed = -sum(pvalue_fixed+log(2), na.rm=TRUE),
                p = min(p),
                p_fixed = min(p_fixed)
            )
            df = as.data.frame(df)
            df <- remove.homozygous.snps(df, blacklist[[ident]], FALSE)
            if (is.null(p_stat)) {
                p_stat <- data.frame(p=cut(df[,'p'], breaks=seq(0, 1, 0.1)), ident=ident, time=time, cov=df[,'mean_cov'])
            } else {
                p_stat <- rbind(p_stat, data.frame(p=cut(df[,'p'], breaks=seq(0, 1, 0.1)), ident=ident, time=time, cov=df[,'mean_cov']))
            }
        }
    }
    plot.p.hist(p_stat, paste0(header, '_all_p_hist'))
}
ranking.pvalue.and.correlation <- function(all_snp, source, header) {
    all_data <- NULL
    for (ident in 1:5) {
        for (time in 1:length(source[[ident]])) {
            print(c(ident, time))
            # print(head(source[[ident]][[time]]))
            df = source[[ident]][[time]] %>% group_by(gene) %>% summarise(
                mean_cov = mean(cov),
                cor      = min(cor),
                scor     = min(scor),
                prob     = -sum(prob, na.rm=TRUE),
                pvalue   = -sum(pvalue+log(2), na.rm=TRUE),
                prob_fixed   = -sum(prob_fixed, na.rm=TRUE),
                pvalue_fixed = -sum(pvalue_fixed+log(2), na.rm=TRUE),
                p = min(p),
                p_fixed = min(p_fixed)
            )
            df = as.data.frame(df)
            df <- remove.homozygous.snps(df, blacklist[[ident]], FALSE)
            for (ase in c(0.5, 0.6, 0.7, 0.8, 0.9, -1)) {
                print(ase)
                tail = paste0('_', ase)
                if (ase < 0) {
                    tail = ''
                    temp <- df[,]
                } else {
                    temp <- df[(df[,'p'] >= ase) & (df[,'p'] < (ase+0.1)),]
                    temp <- rbind(temp, df[(1-df[,'p'] >= ase) & (1-df[,'p'] < (ase+0.1)),])
                }
                print(c(ase, dim(temp)[1]))
                for (y in c('cor', 'scor')[1:1]) {
                    for (x in c('prob', 'pvalue')) {
                        data <- cbind(x=temp[,x], y=temp[,y])
                        data <- data[is.finite(data[,'y']) & is.finite(data[,'x']),]
                        plot.rank.mean.ave(data, paste0(header, '_', ident, '_', time, '_', x, '_', y, tail), all=(ase<0))
                        plot.rank.mean.ave(data, paste0(header, '_', ident, '_', time, '_', x, '_', y, tail, '_rank'), TRUE, all=(ase < 0))
                        colnames(data) <- c('y', 'x')
                        plot.rank.mean.ave(data, paste0(header, '_', ident, '_', time, '_', y, '_', x, tail), all=(ase<0))
                        plot.rank.mean.ave(data, paste0(header, '_', ident, '_', time, '_', y, '_', x, tail, '_rank'), TRUE, all=(ase < 0))
                        colnames(data) <- c('x', 'y')
                        if (y != 'cor' || x != 'pvalue') next
                        if (ase == 'all') {
                            if (is.null(all_data)) all_data <- data
                            else all_data <- rbind(all_data, data)
                        }
                    }
                }
            }
        }
    }
    x <- "pvalue"
    y <- "cor"
    plot.rank.mean.ave(all_data, paste0(header, '_all_', x, '_', y))
    plot.rank.mean.ave(all_data, paste0(header, '_all_', x, '_', y, '_rank'), TRUE)
    colnames(all_data) <- c('y', 'x')
    plot.rank.mean.ave(all_data, paste0(header, '_all_', y, '_', x))
    plot.rank.mean.ave(all_data, paste0(header, '_all_', y, '_', x, '_rank'), TRUE)
}

hist.correlation.grouped <- function(sig_snp, all_snp, source, header) {
    print('correlation')
    print(head(sig_snp))
    print(head(all_snp))
    print(head(source[[1]][[1]]))
    for (ident in 1:5) {
        for (time in 1:length(source[[ident]])) {
            print(c(ident, time))
            df = source[[ident]][[time]] %>% group_by(gene) %>% summarise(
                mean_cov = mean(cov),
                cor      = min(cor),
                scor     = min(scor),
                prob     = sum(prob, na.rm=TRUE),
                pvalue   = sum(pvalue, na.rm=TRUE),
                prob_fixed   = sum(prob_fixed, na.rm=TRUE),
                pvalue_fixed = sum(pvalue_fixed, na.rm=TRUE),
                p = min(p),
                p_fixed = min(p_fixed)
            )
            df = data.frame(df)
            print(head(df))
            for (sig in c('', '_xscaf')[1]) {
                scores <- c('cor', 'scor')
                label = 'sig'
                if (sig != '') label = 'xscaf'
                if (sig == '_xscaf')
                    scores <- c(scores, 'pvalue', 'prob', 'pvalue_fixed', 'prob_fixed')
                scores <- c('pvalue')
                for (cor in scores) { # sig and all
                    print(cor)
                    all_cor <- data.frame(cor=df[,cor], snp=rep('all', dim(df)[1]))
                    if (!any(cor == c('cor', 'scor'))) {
                        all_cor[,"cor"] <- pmax(-50, all_cor[,"cor"])
                    }
                    print(dim(all_cor))
                    if (sig == '') {
                        tsig <- sig_snp[(sig_snp[,"ident"] == ident) & (sig_snp[,"time"] == time),]
                        temp <- subset(df, as.numeric(df[,'gene']) %in% as.numeric(tsig$'gene'))
                        # print(as.numeric(tsig$'gene'))
                    } else {
                        x_scaf <- all_snp[(all_snp[,"ident"] == ident) & (all_snp[,"time"] == time),]
                        x_scaf <- subset(x_scaf, x_scaf$chrm %in% scaffoldsX.prop3)
                        print(head(x_scaf))
                        temp <- subset(df, as.numeric(df[,'gene']) %in% as.numeric(x_scaf$'gene'))
                        print(head(temp))
                        # exit()
                    }
                    if (dim(temp)[1] > 0) 
                        plot.hist.two.correlation(temp, all_cor, cor, paste0(header, '_', ident, '_', time, '_', cor, sig), label)
                    print(c(ident, time, cor))
                    if (!any(cor == c('cor', 'scor')))
                        next
                    if (sig == '') {
                        asig <- sig_snp[(sig_snp[,"ident"] == ident) & (sig_snp[,"time"] <= time),]
                        # tsig <- sig_snp[(sig_snp[,"ident"] == ident) & (sig_snp[,"time"] == time | sig_snp[,"time"] == 1),]
                        tsig <- asig[!duplicated(tsig$'gene'),]
                        tsig <- tsig[unlist(sapply(tsig$'gene', function(x){return(length(which(x == asig$'gene')) == time)})),]
                        print('duplicated')
                        print(dim(tsig))
                        temp <- subset(df, as.numeric(df[,'gene']) %in% as.numeric(tsig$'gene'))
                    } else {
                        x_scaf <- all_snp[(all_snp[,"ident"] == ident) & (all_snp[,"time"] <= time),]
                        x_scaf <- subset(x_scaf, x_scaf$chrm %in% scaffoldsX.prop3)
                        x_scaf <- x_scaf[unlist(sapply(x_scaf$'gene', function(x){return(length(which(x == x_scaf$'gene')) == time)})),]
                        x_scaf <- x_scaf[!duplicated(x_scaf$'gene'),]
                        temp <- subset(df, as.numeric(df[,'gene']) %in% as.numeric(x_scaf$'gene'))
                    }
                    if (dim(temp)[1] > 0) 
                        plot.hist.two.correlation(temp, all_cor, cor, paste0(header, '_', ident, '_', time, '_', cor, sig, '_dup'), label)
                }
            }
            next
            for (cor in c('cor', 'scor')) { # dist with cor-based bins
                for (score in c('prob', 'pvalue')) {
                    all_cor <- data.frame(cor=cut(df[,cor], breaks=20), score=pmax(-50, df[,score]))
                    all_cor <- cbind(all_cor, pn=unlist(sapply(df[,cor], function(x){if(is.na(x)) return(NA); if(x > 0) return('positive'); if (x < 0) return('negative'); return('zero')})))
                    all_cor <- all_cor[!is.na(all_cor[,'cor']),]
                    plot.boxplot.correlation(all_cor, paste0(cor, '_', score, '_', ident, '_', time))
                }
            }
        }
    }
}

compute.mean.value <- function(df, ident) {
    tdf =  df %>% group_by(gene, ensemblID, chrm, pos) %>% summarise(
        mean_cov = mean(mean_cov),
        mean_pvalue = mean(pvalue, na.rm=TRUE),
        max_pvalue  = min(pvalue, na.rm=TRUE),
        sd_pvalue   = sd(pvalue, na.rm=TRUE),
        mean_prob = mean(prob, na.rm=TRUE),
        max_prob  = min(prob, na.rm=TRUE),
        sd_prob   = sd(prob, na.rm=TRUE),
        count = n()
    )
    tdf <- cbind(tdf, ident=rep(ident, dim(tdf)[1]))
    return(tdf)
}

filter.and.factor.conversion <- function(df, thres=50) {
    df$count = factor(df$count, levels=c(1, 2, 3))
    df$mean_pvalue = pmax(-thres, df$mean_pvalue)
    df$max_pvalue = pmax(-thres, df$max_pvalue)
    df$mean_prob = pmax(-thres, df$mean_prob)
    df$max_prob = pmax(-thres, df$max_prob)
    return(df)
}

compute.mean.ident.value <- function(df) {
    tdf <- df %>% group_by(ensemblID, chrm, pos) %>% summarise(
        all_count = sum(count),
        count = n(),
        mean_mean_cov = mean(mean_cov),
        median_mean_pvalue = median(mean_pvalue, na.rm=TRUE),
        max_mean_pvalue = min(mean_pvalue, na.rm=TRUE),
        sd_mean_pvalue = sd(mean_pvalue, na.rm=TRUE),
        median_mean_prob = median(mean_prob, na.rm=TRUE),
        max_mean_prob = min(mean_prob, na.rm=TRUE),
        sd_mean_prob = sd(mean_prob, na.rm=TRUE),
    )
    return(tdf)
}
filter.and.factor.conversion.for.mean <- function(df, thres=50) {
    df$count = factor(df$count)
    df$median_mean_pvalue = pmax(-50, df$median_mean_pvalue)
    df$max_mean_pvalue = pmax(-50, df$max_mean_pvalue)
    df$median_mean_prob = pmax(-50, df$median_mean_prob)
    df$max_mean_prob = pmax(-50, df$max_mean_prob)
    return(df)
}

compute.recurrent.significance <- function(all_snp, source, header) {
    require(ggplot2)
    require(ggsci)
    require(dplyr)
    print(head(all_snp))
    print(head(source[[ident]][[1]]))
    return()
    # all_df = NULL
    # for (ident in 1:5) {
    #     temp = all_snp[all_snp[,"ident"] == ident,]
    #     print(head(temp))
    #     temp[,"pos"] <- as.character(temp[,"pos"])
    #     for (time in 1:3) {

    #     }
    #     df =  compute.mean.value(temp, ident)
    #     stopifnot(max(df$count) <= 3)
    #     if (is.null(all_df)) all_df <- df
    #     else all_df <- rbind(all_df, df)
    #     df = filter.and.factor.conversion(df)
    #     for (y in c('max', 'sd')) {
    #         for (score in c('pvalue', 'prob')) {
    #             png(paste0('mean_', y, '_each_ind_', ident, '_', header, '_', score, '.png'))
    #             g <- ggplot(data=df, aes_string(x=paste0('mean_', score), y=paste0(y, '_', score)))+geom_point(aes(color=count), alpha=0.3)+theme_bw()+scale_color_manual(values=c("1"="gray", "2"="blue", "3"="red"))
    #             plot(g)
    #             dev.off()
    #         }
    #     }
    # }
    # print(head(all_df))
    # df = compute.mean.ident.value(all_df)
    # df = filter.and.factor.conversion(df)
    # for (y in c('max', 'sd')) {
    #     for (score in c('pvalue', 'prob')) {
    #         png(paste0('median_', y, '_all_ind_', header, '_', score, '.png'))
    #         g <- ggplot(data=df, aes_string(x=paste0('median_mean_', score), y=paste0(y, '_mean_', score)))+geom_point(aes(color=count), alpha=0.5)+theme_bw()+scale_color_npg()
    #         plot(g)
    #         dev.off()
    #     }
    # }
}
plot.pvalues.grouped <- function(all_snp, source, header) {
    require(ggplot2)
    require(ggsci)
    require(dplyr)
    print(head(all_snp))
    all_df = NULL
    for (ident in 1:5) {
        temp = all_snp[all_snp[,"ident"] == ident,]
        print(head(temp))
        temp[,"pos"] <- as.character(temp[,"pos"])
        df =  compute.mean.value(temp, ident)
        stopifnot(max(df$count) <= 3)
        if (is.null(all_df)) all_df <- df
        else all_df <- rbind(all_df, df)
        df = filter.and.factor.conversion(df)
        for (y in c('max', 'sd')) {
            for (score in c('pvalue', 'prob')) {
                png(paste0('mean_', y, '_each_ind_', ident, '_', header, '_', score, '.png'))
                g <- ggplot(data=df, aes_string(x=paste0('mean_', score), y=paste0(y, '_', score)))+geom_point(aes(color=count), alpha=0.3)+theme_bw()+scale_color_manual(values=c("1"="gray", "2"="blue", "3"="red"))
                plot(g)
                dev.off()
            }
        }
    }
    print(head(all_df))
    df = compute.mean.ident.value(all_df)
    df = filter.and.factor.conversion(df)
    for (y in c('max', 'sd')) {
        for (score in c('pvalue', 'prob')) {
            png(paste0('median_', y, '_all_ind_', header, '_', score, '.png'))
            g <- ggplot(data=df, aes_string(x=paste0('median_mean_', score), y=paste0(y, '_mean_', score)))+geom_point(aes(color=count), alpha=0.5)+theme_bw()+scale_color_npg()
            plot(g)
            dev.off()
        }
    }
}

overlap.sig.snps <- function(sig_snp, all_snp, all_ratio) {
    print('correlation')
    for (ident in 1:5) {
        sig_snp_genes <- lapply(0:2, function(x) {
            return(sig_snp[(sig_snp[,"ident"] == ident) & (sig_snp[,"time"] == x+1), 'gene'])
        })
        for (time in 2:length(sig_snp_genes)) {
            all_genes = subset(all_snp, (all_snp[,"ident"] == ident) & (all_snp[,"time"] == 1 | all_snp[,"time"] == time))[,"gene"]
            all_genes = length(unique(all_genes))
            int_genes = intersect(sig_snp_genes[[1]], sig_snp_genes[[time]])
            N = length(sig_snp_genes[[1]])
            M = length(sig_snp_genes[[time]])
            n = length(int_genes)
            # print(c(n, N, M, all_genes))
            data <- matrix(c(n, N-n, M-n, all_genes-M-N+n), nrow=2)
            # pvalue <- hypergeom.test(data, alternative='greater')
            # print(phyper(n-1, M, all_genes-M, N, lower.tail=TRUE, log.p=FALSE))
            pvalue <- 1.0-phyper(n-1, M, all_genes-M, N, lower.tail=TRUE, log.p=FALSE)
            print(c(n, N, M, all_genes, pvalue))
            if (time == 2 && time != length(sig_snp_genes)) {
                int_genes = intersect(intersect(sig_snp_genes[[1]], sig_snp_genes[[time]]), sig_snp_genes[[time+1]])
                print(c(ident, 'triple', length(int_genes), N, M, all_genes))
            }
        }
    }
}

check.raw.ase.dist <- function(rdf, all_df, header) {
    require(ggplot2)
    require(viridis)
    require(ggsci)
    require(dplyr)
    for (prop in c(1:3)) {
        scaffold = scaffoldsX.prop
        if (prop == 2) scaffold = scaffoldsX.prop2
        if (prop == 3) scaffold = scaffoldsX.prop3
        for (time in 1:3) {
            for (ident in 1:5) {

                print(head(rdf[[ident]][[time]]))
                print(dim(rdf[[ident]][[time]]))
                temp = all_df[(all_df[,"ident"] == ident) & (all_df[,"time"] == time),]
                m <- merge(rdf[[ident]][[time]][,c('x', 'gene')], temp[,c("gene", "chrm")], by="gene")
                # m <- m[1:10000,]
                print(head(m))
                m[,"x"] = pmax(m[,"x"], 1-m[,"x"])
                # m <- m %>% group_by(gene, chrm)  %>% summarise(
                #     x = mean(x, na.rm=TRUE)
                # )
                # m <- as.data.frame(m)
                print(head(m))
                print(dim(m))
                png(paste0('hist_ase_', ident, '_', time, '_all_', header, '_', prop, '.png'))
                g <- ggplot(data=m, aes(x=x, group=chrm))+geom_freqpoly(color='gray', alpha=0.5)+theme_bw()
                plot(g)
                dev.off()
                chrom <- unique(m[,"chrm"])
                m[,"type"] <- unlist(sapply(m[,"chrm"], function(x){ if(any(x == scaffold))return('red'); return('gray')}))
                m <- m[order(m[,"type"]),]
                m[,"chrm"] <- factor(m[,"chrm"])
                print(head(m))
                col_dict <- sapply(unique(m[,"chrm"]), function(x){if(any(as.character(x) == scaffold))return('red'); return('gray')})
                print(head(col_dict))
                png(paste0('hist_ase_', ident, '_', time, '_all_xcol_', header, '_', prop, '.png'))
                g <- ggplot(data=m, aes(x=x, color=chrm))+geom_freqpoly()+theme(legend.position="none")+scale_color_manual(values=col_dict)
                plot(g)
                dev.off()
            }
        }
    }
}

args = commandArgs(trailingOnly=TRUE)
print(args)
if (args[1] == "hist")
    hist.correlation.grouped(sig_df, all_df, sdf, 'each_null')

# plot.pvalues.grouped(all_df, sdf, 'each_null')
# compute.recurrent.significance(all_df, sdf, 'each_null')
# check.raw.ase.dist(rdf, all_df, 'each_null')
if (args[1] == "rank")
    ranking.pvalue.and.correlation(all_df, sdf, 'each_null')
if (args[1] == 'phist')
    plot.basic.statistics(all_df, rdf, 'each_null')
# overlap.sig.snps(sig_df, all_df)
# hist.correlation.grouped(sig_t0, all_t0, sdf, 't0_null')
# binned.correlation(sdf, sdf, 't0_null')
# binned.correlation(rdf, sdf, 'each_null')

