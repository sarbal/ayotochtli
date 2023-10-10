require(ggplot2)
require(ggpubr)
require(reshape2)
require(viridis)
require(parallel)
sdf <- readRDS('simple_t0_result.rds')

read.fisher.data <- function(fix, i, remx, label) {
    fname = paste0("fisher_stat_obs_emp", fix, i, remx, '.tsv')
    print(fname)
    if (!file.exists(fname)) 
        return(NULL)
    a <- read.table(fname, header=T, sep="\t")
    print(head(a))
    ident = floor((i-1)/3+1)
    time = ((i-1)%%3+1)
    print(dim(a))
    if (label == 'p')
        size <- dim(a)[1]
    else {
        if (any(label == c('emp_p_adj', 'chi_p_adj')))
            a[,label] <- -log10(a[,label])
        size <- length(which(a[,label] > -log10(0.05)))
    }
    if (fix != '_') {
        if (label == 'p')
            return(list(size, a[,'p_fixed'], a[,'gene']))
        else
            return(list(size, a[,label], a[,'gene']))
    } else {
        return(list(size, a[,label], a[,'gene']))
    }
}

plot.hist.all.snp.groups <- function() {
    for (cor in c('cor', 'scor')) {
        for (fix in c('_', '_fixed_')) {
            for (label in c('emp_p_adj')) {
            # for (label in c('p', 'emp_p', 'chi_p', 'emp_p_adj', 'chi_p_adj')) {
                all_table <- NULL
                for (remx in c('', '_remx', '_extx')) {
                    size = c()
                    pdata = NULL
                    for (i in 1:15) { 
                        result <- read.fisher.data(fix, i, remx, label)
                        if (is.null(result)) next
                        size <- c(size, result[[1]])
                        ident = floor((i-1)/3+1)
                        time = ((i-1)%%3+1)
                        rtemp <- cbind(value=result[[2]], gene=result[[3]])
                        if ((i-1) %% 3 == 0) {
                            rtemp <- cbind(rtemp, rep(NA, dim(rtemp)[1]))
                            colnames(rtemp) <- c('value', 'gene', 'cor')
                        } else {
                            temp = sdf[[floor((i-1)/3+1)]][[(i-1)%%3]]
                            temp <- temp[,c('gene', cor)]
                            temp <- temp[!duplicated(temp),]
                            colnames(temp) <- c('gene', 'cor')
                            rtemp <- merge(rtemp, temp, by='gene')
                        }
                        print(head(rtemp))
                        if (is.null(pdata)) pdata <- rtemp
                        else pdata <- rbind(pdata, rtemp)
                    }
                    pdata <- cbind(pdata, snp=rep(remx, dim(pdata)[1]))
                    if (is.null(all_table)) all_table <- pdata
                    else all_table <- rbind(all_table, pdata)
                }
                all_table <- all_table[!is.na(all_table[,'cor']),] # remove first time point                           
                all_table[,4] <- sapply(all_table[,4], function(x){if(x == "") return("All"); if(x == "_remx") return("No X"); return("X");})
                png(paste0('violinplot_ase_', fix, '_', '_', cor, '.png'))
                plot(ggviolin(all_table, x='snp', y='cor', color='snp', add = "boxplot"))
                dev.off()
                df <- all_table %>%   
                        group_by(snp) %>% 
                        mutate(rank = rank(value, ties.method = "first"))
                df <- data.frame(df)
                count <- table(df[,'snp'])
                df[,'rank'] <- unlist(mclapply(1:dim(df)[1], function(x){return(df[x,'rank']/count[names(count) == df[x,'snp']])}, mc.cores=50))
                write.table(df, paste0('rank_data_', fix, '_all_', label, '_', cor, '.csv'), sep=",")
                methods <- c('glm', 'gam', 'lm')
                colnames(df)[4] = "SNP"
                for (method in methods) {
                    ylab = "Pearson Correlation Coefficient"
                    if (cor == 'scor') ylab = "Spearman Correlation Coefficient"
                    xlab = "Ranking order of -log10 FDR"
                    p <- ggplot(df, aes(x=rank, y=cor, color=SNP) ) +
                            geom_smooth(aes(color=SNP), method=method) + theme_bw()+theme(text = element_text(size=20)) + xlab(xlab)+ylab(ylab)
                    pdf(paste0('rankplot_binned_p_', fix, '_all', '_', label, '_', cor, '_', method, '.pdf'))
                    plot(p)
                    dev.off()
                }
            }
        }
    }
}

plot.hist.each.snp.groups <- function() {
    for (cor in c('cor', 'scor')) {
        for (fix in c('_', '_fixed_')) {
            for (remx in c('', '_remx', '_extx')) {
                for (label in c('p', 'emp_p', 'chi_p', 'emp_p_adj', 'chi_p_adj')) {
                    pdata = NULL
                    size = c()
                    if (FALSE) {
                        for (i in 1:15) { 
                            ident = floor((i-1)/3+1)
                            time = ((i-1)%%3+1)
                            result <- read.fisher.data(fix, i, remx, label)
                            if (is.null(result)) next
                            size <- c(size, result[[1]])
                            pdata <- cbind(pdata, result[[2]])
                            colnames(pdata)[dim(pdata)[2]] <- paste0(ident, '_', time)
                            names(size)[length(size)] <- paste0(ident, '_', time)
                        }
                        print(head(pdata))
                        # colnames(pdata) = paste0((1:15-1)/5+1, '_', (1:15-1)%%3+1)
                        pdata <- melt(pdata)[,c(2,3)]
                        colnames(pdata) <- c('label', 'value')
                        png(paste0('plot_ase_', fix, remx, '_', label, '.png'))
                        if (label == 'p')
                            plot(ggdensity(pdata, x = 'value', rug=FALSE, color="label", xlim=c(0, 1)))
                        else
                            plot(ggdensity(pdata, x = 'value', rug=FALSE, color="label"))
                        dev.off()
                        size <- data.frame(value=size, label=names(size))
                        print(head(size))
                        png(paste0('snp_number_ase_', fix, remx, '_', label, '.png'))
                        plot(gghistogram(size, x = 'value', color="label", fill="label"))
                        dev.off()
                        print(head(pdata))
                    } else {
                        pdata <- NULL
                        for (i in 1:15) { 
                            result <- read.fisher.data(fix, i, remx, label)
                            if (is.null(result)) next
                            rtemp <- cbind(value=result[[2]], gene=result[[3]])
                            if ((i-1) %% 3 == 0) {
                                rtemp <- cbind(rtemp, rep(NA, dim(rtemp)[1]))
                                colnames(rtemp) <- c('value', 'gene', 'cor')
                            } else {
                                print(c(floor((i-1)/3+1), (i-1)%%3))
                                temp = sdf[[floor((i-1)/3+1)]][[(i-1)%%3]]
                                temp <- temp[,c('gene', cor)]
                                temp <- temp[!duplicated(temp),]
                                colnames(temp) <- c('gene', 'cor')
                                rtemp <- merge(rtemp, temp, by='gene')
                                print('merged')
                                print(head(rtemp))
                            }
                            if (is.null(pdata)) pdata <- rtemp
                            else pdata <- rbind(pdata, rtemp)
                        }
                        pdata <- pdata[!is.na(pdata[,'cor']),]
                        p <- pdata %>%               
                            mutate(bin=cut_width(value, width=max(pdata[,'value'], na.rm=TRUE)/15, boundary=0) ) %>%
                            ggplot( aes(x=bin, y=cor) ) +
                                geom_boxplot(fill='gray') + theme_bw() + xlab(paste0("-log10 ", label)) + geom_abline(intercept = 0, slope = 0, color='darkgray')
                        png(paste0('boxplot_binned_p_', fix, remx, '_', label, '_', cor, '.png'))
                        plot(p)
                        dev.off()
                        p <- pdata %>%               
                            mutate(rank=rank(pdata[,'value'])) %>%
                            ggplot( aes(x=rank, y=cor) ) +
                                geom_point(alpha=0.1) + geom_smooth() + theme_bw() + xlab(paste0("-log10 ", label, ' rank order'))
                        png(paste0('rankplot_binned_p_', fix, remx, '_', label, '_', cor, '.png'))
                        plot(p)
                        dev.off()
                        png(paste0('scatter_p_', fix, remx, '_', label, '_', cor, '.png'))
                        p <- ggplot(pdata, aes(x=value, y=cor)) + geom_point(alpha=0.05)+theme_bw() + xlab(paste0("-log10 ", label))
                        plot(p)
                        dev.off()
                    }
                }
            }
        }
    }
}



plot.hist.all.snp.groups()
# plot.hist.each.snp.groups()