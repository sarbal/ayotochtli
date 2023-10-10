# simulate the binomial accuracy and coverage and cell stage
require(ggpubr)
require(dplyr)
require(tidyr)
require(viridis)
require(heatmaply)
require(reshape2)
require(ComplexHeatmap)
require(matrixStats)
require(ggsci)
require(cowplot)
require(pROC)
require(ramify)
require(parallel)

chi.test <- function(ase, cov) {
    ma <- which.max(ase)
    mi <- which.min(ase)
    tbl <- rbind(c(ase[ma], cov-ase[ma]), c(ase[mi], cov-ase[mi]))
    if (min(tbl) <= 5) { # Fisher
        return(fisher.test(tbl)$p.value)
    } else {
        return(chisq.test(tbl)$p.value)
    }
}
comp.zscore.p <- function(ase, index) {
    z <- mean(ase[-index])-ase[index]/(sd(ase[-index])/sqrt(3))
    return(2*pnorm(-abs(z)))
}
zscore.test <- function(ase) {
    ma <- which.max(ase)
    mi <- which.min(ase)
    p.value <- min(comp.zscore.p(ase, ma), comp.zscore.p(ase, mi))
    return(p.value)
}

plot.smooth.roc.data <- function(header, data) {
    print('roc')
    for (i in c(10, 25)) {
        df <- data
        df[,'FPR'] <- unlist(sapply(df[,'FPR'], function(x){return(ceiling(x*i)/i)}))
        # df[,'FPR'] <- cut(df[,'FPR'], i)
        df <- df %>% group_by(quad, FPR) %>% summarise(
            sd_TPR=sd(TPR),
            TPR=mean(TPR),
        )
        df <- as.data.frame(df)
        quad <- unique(df[,'quad'])
        df <- rbind(data.frame(quad=quad, FPR=rep(0, length(quad)), sd_TPR=rep(1, length(quad)), TPR=rep(0, length(quad))), df)
        df[,'FPR'] <- factor(df[,'FPR'])
        ymax =max(1, max(df[,'TPR']+df[,'sd_TPR']))
        # labels <- paste0('-', unlist(gsub(']', '', unlist(sapply(unique(df[,'FPR']), function(x){return(strsplit(as.character(x), ',')[[1]][[2]])})))))
        labels <- as.numeric(as.character(unique(df[,'FPR'])))
        names(labels) <- unique(df[,'FPR'])
        print(labels)
        g <- ggplot(df, aes(x=FPR, y=TPR, color=quad))+theme_cowplot()+xlab('FPR')+ylab('TPR')
        g <- g+geom_point(position=position_dodge(.9))+ geom_linerange(aes(ymin=TPR-sd_TPR, ymax=TPR+sd_TPR), position=position_dodge(.9))
        g <- g+scale_x_discrete('FPR', labels=labels)+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+scale_color_npg()+coord_cartesian(ylim=c(0, ymax))
        pdf(paste0(header, '_msd_', i, '.pdf'))
        plot(g)
        dev.off()
        g <- ggplot(df, aes(x=FPR, y=TPR, color=quad))+theme_cowplot()+xlab('FPR')+ylab('TPR')
        g <- g+geom_point(size=3)+geom_path(aes(group=quad), lwd=1.5)
        g <- g+scale_x_discrete('FPR', labels=labels)+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+scale_color_npg()+coord_cartesian(ylim=c(0, 1))
        pdf(paste0(header, '_msd_', i, '_mean.pdf'))
        plot(g)
        dev.off()
    }
    for (i in c(3)) {
        var <- as.formula(paste0('y ~ poly(x,', i, ')'))
        g <- ggplot(data, aes(x=FPR, y=TPR, color=quad))+theme_cowplot()+xlab('FPR')+ylab('TPR')
        g <- g +  geom_smooth(method = 'lm', formula = var, span=0.7, se = TRUE)+scale_color_npg()+coord_cartesian(ylim=c(0, ymax))
        pdf(paste0(header, '_smo_', i, '.pdf'))
        plot(g)
        dev.off()
    }
}

compute.auroc <- function(cov, result, SNP, ratio, times) {
    all_data <- NULL
    each_roc <- NULL
    all_auc <- NULL
    y_vector <- c(rep(1, as.integer(SNP*ratio)), rep(0, as.integer(SNP-SNP*ratio)))
    for (i in seq(1, 20)) {
        p <- matrix(p.adjust(as.vector(result[1:i,]), method='BH'), ncol=SNP)
        p <- colMins(p)
        roc_obj <- roc(y_vector, -p)
        roc_df <- data.frame(
            TPR=rev(roc_obj$sensitivities), 
            FPR=rev(1 - roc_obj$specificities)
        )
        temp <- cbind(cbind(roc_df$FPR, roc_df$TPR), rep(i, dim(roc_df)[1]))
        if (is.null(each_roc)) each_roc <- temp
        else each_roc <- rbind(each_roc, temp)
        all_auc <- rbind(all_auc, c(cell, cov, i, auc(roc_obj)))
        TP <- length(which(p[1:as.integer(SNP*ratio)] <= 0.05))
        FP <- length(which(p[as.integer(SNP*ratio+1):SNP] <= 0.05))
        print(c(TP, FP, as.integer(SNP*ratio), as.integer(SNP-SNP*ratio)))
        stopifnot(TP <= SNP*ratio && FP <= SNP-SNP*ratio)
        TP <- min(TP/as.integer(SNP*ratio), 1.0)
        FP <- min(FP/as.integer(SNP-SNP*ratio), 1.0)
        stopifnot(TP <= 1 && FP >= 0)
        all_data <- rbind(all_data, c(i, TP, FP, auc(roc_obj), times, cov))
    }
    return(list(all_data, each_roc, all_auc))
}

plot.smooth.tp <- function(all_data, header) {
    print('tp')
    ylab = "True Positive rate"
    all_data[,'FPR'] <- clip(all_data[,'FPR'], 0., 1.)
    # xlab = "1-False Positive rate"
    for (xax in c('FPR', 'quad')) {
        if (xax == 'quad') xlab = "Quad"
        else xlab = "False Positive rate"
        for (i in c(10, 25)) {
            if (xax == 'quad') {
                df <- all_data
                df[,'quad'] <- factor(df[,'quad'])
                df <- df %>% group_by(quad, cov) %>% summarise(
                    sd_TPR=sd(TPR),
                    TPR=mean(TPR),
                )
                df <- as.data.frame(df)
                labels <- unique(df[,'quad'])
                names(labels) <- unique(df[,'quad'])
            } else {
                df <- all_data
                df[,'FPR'] <- cut(df[,'FPR'], i)
                df <- df %>% group_by(FPR, cov) %>% summarise(
                    sd_TPR=sd(TPR),
                    TPR=mean(TPR),
                )
                df <- as.data.frame(df)
                labels <- paste0('-', unlist(gsub(']', '', unlist(sapply(unique(df[,'FPR']), function(x){return(strsplit(as.character(x), ',')[[1]][[2]])})))))
                names(labels) <- unique(df[,'FPR'])
            }
            ymax =max(1, max(df[,'TPR']+df[,'sd_TPR']))+0.08
            g <- ggplot(df, aes_string(x=xax, y='TPR', color='cov'))+theme_cowplot()+xlab(xlab)+ylab(ylab)
            g <- g+geom_point(position=position_dodge(.9))+ geom_linerange(aes(ymin=TPR-sd_TPR, ymax=TPR+sd_TPR), position=position_dodge(.9))+ coord_cartesian(ylim=c(0, ymax))
            g <- g+scale_x_discrete(xlab, labels=labels)+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+scale_color_npg()
            pdf(paste0('plot_', xax, '_', header, '_msd_', i, '.pdf'))
            plot(g)
            dev.off()
        }
        for (i in c(3)) {
            var <- as.formula(paste0('y ~ poly(x,', i, ')'))
            g <- ggplot(all_data, aes_string(x=xax, y='TPR', color='cov'))+theme_cowplot()+xlab(xlab)+ylab('TPR')
            g <- g +  geom_smooth(method = 'lm', formula = var, span=0.7, se = TRUE)+scale_color_npg()
            pdf(paste0('plot_', xax, '_', header, '_smo_', i, '.pdf'))
            plot(g)
            dev.off()
        }
        print('end')
    }
}
simulate.pvalue.auc <- function(SNP, ratio, mtimes, cov, cov_set, method, cell, quad) 
{
    all_auc <- NULL
    all_data <- NULL
    cheader = paste0(cell, '_', mtimes, '_', ratio, '_', method, '_', cov)
    each_roc <- NULL
    for (times in 1:mtimes) {
        prob_samples <- rbinom(quad*4, cell, 0.5)/cell
        ase_vec = prob_samples
        result <- do.call(rbind, lapply(1:quad, function(x) {
            start = (x-1)*4+1
            end = (x-1)*4+4
            p <- c()
            for (snp in 1:as.integer(SNP*ratio)) {
                ase = sapply(prob_samples[start:end], function(x){return(rbinom(1, cov, x))})
                p <- c(p, chi.test(ase, cov))
            }
            return(p)
        }))
        cont_result <- do.call(rbind, lapply(1:quad, function(x) {
            return(runif((SNP-as.integer(SNP*ratio)), min=0, max=1))
        }))
        result <- as.matrix(cbind(result, cont_result))
        roc_result <- compute.auroc(cov, result, SNP, ratio, times)
        all_data <- rbind(all_data, roc_result[[1]])
        each_roc <- rbind(each_roc, roc_result[[2]])
        all_auc <- rbind(all_auc, roc_result[[3]])
    }
    colnames(each_roc) <- c('FPR', 'TPR', 'quad')
    each_roc <- as.data.frame(each_roc)
    stopifnot(max(each_roc[,'TPR']) <= 1.0)
    saveRDS(each_roc, paste0('each_roc_', cheader, '.rds'))
    each_roc <- each_roc[sapply(each_roc[,'quad'], function(x){return(any(x == c(1, 2, 5, 10, 20)))}), ]
    each_roc[,'quad'] <- factor(each_roc[,'quad'])
    if (any(cov == cov_set))
        plot.smooth.roc.data(paste0('plot_roc_', cheader), each_roc)
    return(list(all_data, all_auc))
}

plot.auroc.only <- function(SNP, ratio, mtimes, cov, cov_set, method, cell, max_quad) {
    cheader = paste0(cell, '_', mtimes, '_', ratio, '_', method, '_', cov)
    each_roc <- readRDS(paste0('each_roc_', cheader, '.rds'))
    each_roc <- each_roc[sapply(each_roc[,'quad'], function(x){return(any(x == c(1, 2, 5, 10, 20)))}), ]
    each_roc[,'quad'] <- factor(each_roc[,'quad'])
    print(cheader)
    plot.smooth.roc.data(paste0('plot_roc_', cheader), each_roc)
}

max_quad = 100
set.seed(4)
SNP <- 10000
mtimes <- 100
CORES <- 20
method <- "lm"
init <- 10
cov_set <- c(10, 25, 50, 75, 100)
if (FALSE) {
    for (ratio in c(0.001, 0.01, 0.05)) {
        all_auc <- NULL
        aheader = paste0(mtimes, '_', ratio, '_', method)
        for (cell in c(2, 4, 8, 16, 32, 64, 128, 256, 512, 10000)) {
            all_data <- NULL
            header = paste0(cell, '_', mtimes, '_', ratio, '_', method)
            all_data <- NULL
            results <- mclapply(init:200, function(cov) {return(simulate.pvalue.auc(SNP, ratio, mtimes, cov, cov_set, method, cell, max_quad))}, mc.cores=CORES)
            print(length(results))
            for (i in init:200) {
                all_auc <- rbind(all_auc, results[[i-init+1]][[2]])
                all_data <- rbind(all_data, results[[i-init+1]][[1]])
            }
            colnames(all_data) <- c('quad', 'TPR', 'FPR', 'auc', 'rep', 'cov')
            all_data <- as.data.frame(all_data)
            write.table(all_data, paste0('output_', header, '.txt'))
            all_data <- all_data[unlist(sapply(all_data[,'cov'], function(x){return(x %in% cov_set)})), ]
            all_data[,'cov'] <- as.factor(all_data[,'cov'])
            plot.smooth.tp(all_data, header)
        }
        write.table(all_auc, paste0('output_auc_', aheader, '.txt'))
    }
} else {
    for (ratio in c(0.001, 0.01, 0.05)) {
        aheader = paste0(mtimes, '_', ratio, '_', method)
        for (cell in c(2, 4, 8, 16, 32, 64, 128, 256, 512, 10000)) {
        # for (cell in c(64)) {
            header = paste0(cell, '_', mtimes, '_', ratio, '_', method)
            for (cov in c(20, 50, 75)) {
                plot.auroc.only(SNP, ratio, mtimes, cov, cov_set, method, cell, max_quad)
            }
        }
    }
}
        




