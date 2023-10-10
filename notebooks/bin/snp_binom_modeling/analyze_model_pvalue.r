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
load("data/Xscaffolds.Rdata")

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

# print(head(rdf[[1]][[1]]))
# print(head(sdf[[1]][[1]]))
# print(head(sig_df))
# print(head(sig_t0))
# print(head(all_df))

extract.homozygous.snp <- function() {
    black_list_all_quads <- lapply(1:5, function(x) {
        black_list_overlap <- c()
        for (i in 1:3) {
            homo <- rdf[[x]][[1]] %>% group_by(gene) %>% summarize(
                min_x = min(x),
                max_x = max(x)
            )
            homo <- as.data.frame(homo)
            blacklist <- subset(homo, homo$min_x == homo$max_x & (homo$min_x == 0 | homo$min_x == 1))[,'gene']
            if (i == 1) black_list_overlap <- blacklist
            else {
                black_list_overlap <- c(black_list_overlap, blacklist)
                black_list_overlap <- unique(black_list_overlap)
                # black_list_overlap <- black_list_overlap[duplicated(black_list_overlap)]
            }
        }
        return(black_list_overlap)
    })
    saveRDS(black_list_all_quads, 'black_list_homo_snps.rds')
}

add.rep.cor.label <- function(tdf, label=FALSE) {
    vec <- sapply(1:dim(tdf)[1], function(i) {
            x <- tdf[i, 'first_pair']
            y <- tdf[i, 'second_pair']
            if (is.na(x) || is.na(y)) return('nan')
            if (x > 0 && y > 0) {
                return('pos')
            } else if (x < 0 && y < 0) {
                return('neg')
            } else {
                return('inc')
            }
    })
    if (label) vec <- convert.pn.label(vec)
    return(vec)
}


store.sig.cor <- function(gene_list, sgene, ident, source, header) {
    require(ggplot2)
    require(ggsci)
    require(ggpubr)
    print(ident)
    require(reshape2)
    for (cor_label in c('scor', 'cor')) {
        print(head(source))
        tdf <- source[,c('gene', cor_label, 'time')]
        tdf <- tdf[!duplicated(tdf[,c('gene', 'time')]),]
        print(dim(tdf))
        tdf <- cbind(tdf, sapply(tdf[,'gene'], function(x){
            vec <- c()
            for (n in 1:length(gene_list)){
                if (any(x == gene_list[[n]])) return(n);
            }
            if (length(vec) > 0) return(length(vec))
            return('others')
        }))
        for (i in 1:3) {
            sig_snp_time = subset(sgene, sgene[,'time'] == i-1)
            tdf <- cbind(tdf, unlist(sapply(tdf[,'gene'], function(x) {
                return(any(x == sig_snp_time[,1]))
            })))
        }
        colnames(tdf) <- c('gene', 'cor', 'time', 'significance', 'sig_1', 'sig_2', 'sig_3')
        tdf <- dcast(tdf, gene + significance + sig_1 + sig_2 + sig_3 ~ time, value.var=c('cor'))
        print(head(tdf))
        colnames(tdf) <- c('gene', 'significance', 'sig_1', 'sig_2', 'sig_3', 'first_pair', 'second_pair')
        tdf <- cbind(tdf, alpha=sapply(tdf[,'significance'], function(x){if(x == 'others') return(0.); return(1.0)}))
        tdf <- cbind(tdf, rep=add.rep.cor.label(tdf))
        print(summary(tdf[,'rep']))
        tdf[,'gene'] <- as.numeric(tdf[,'gene'])
        all_df[,'gene'] <- as.numeric(all_df[,'gene'])
        mdf <- merge(tdf, all_df[all_df[,'ident'] == ident & all_df[,'time'] == 1, c('gene', 'ensemblID', 'type', 'name', 'pos')], by='gene', all.x=TRUE)
        stopifnot(!is.na(mdf[1,'ensemblID']))
        write.table(mdf, file=paste0(header, '_', cor_label, '_comp.csv'), sep=',')
    }
    # detach(plyr)
}

sig.snps.across.time <- function(gene_list, sgene, emp_p_df, pvalue, header, ident, threshold=0.05) {
    print(head(sdf[[ident]][[1]]))
    source = cbind(sdf[[ident]][[1]], time=2)
    source = rbind(source, cbind(sdf[[ident]][[2]], time=3))
    source[,'time'] = factor(source[,'time'])
    store.sig.cor(gene_list, sgene, ident, source, header)
}

plot.cov.and.ase <- function(gene_list, ident, header) {
    require(ggplot2)
    require(ggsci)
    require(ggpubr)
    print(ident)
    # print(head(rdf))
    source = cbind(rdf[[ident]][[1]], time=1)
    source = rbind(source, cbind(rdf[[ident]][[2]], time=2))
    source = rbind(source, cbind(rdf[[ident]][[3]], time=3))
    source[,'time'] = factor(source[,'time'])
    print(head(source))
    for (i in 1:length(gene_list)) {
        if (length(gene_list[[i]]) == 0) next
        temp <- source[source[,'gene'] %in% gene_list[[i]],]
        print(head(source))
        for (gene in unique(temp[,'gene'])) {
            for (time in 1:3) {
                print(gene)
                # print(temp[temp[,'gene'] == gene,])
                stopifnot(length(which(temp[,'gene'] == gene & temp[,'time'] == time)) <= 4)
            }
        }
        print(head(temp))
        temp[,'log10_cov'] <- log10(temp[,'cov'])
        p <- ggboxplot(temp, x = "gene", y = "x",
                color = "time", palette =c("#00AFBB", "#E7B800", "#FC4E07"),
                add = "jitter", shape = "time", outlier.shape = NA)+theme(axis.text.x = element_text(angle = 90, hjust = 1))

        pdf(paste0(header, i, '_ase.pdf'), width=16)
        # png(paste0(header, i, '_ase.png'), width=1000)
        plot(p)
        dev.off()
        p <- ggboxplot(temp, x = "gene", y = "log10_cov",
                color = "time", palette =c("#00AFBB", "#E7B800", "#FC4E07"),
                add = "jitter", shape = "time", outlier.shape = NA)+theme(axis.text.x = element_text(angle = 90, hjust = 1))
        pdf(paste0(header, i, '_cov.pdf'), width=16)
        # png(paste0(header, i, '_cov.png'), width=1000)
        plot(p)
        dev.off()
        print(head(temp))
        # exit()
        temp = cbind(temp, rowname=1:dim(temp)[1])
        mdf <- merge(temp, all_df[all_df[,'ident'] == ident & all_df[,'time'] == 1, c('gene', 'ensemblID', 'type', 'name', 'pos')], by='gene', all.x=TRUE)
        mdf <- mdf[order(mdf[,'gene'], as.numeric(mdf[,'rowname'])),]
        mdf <- mdf[, -c(which(colnames(mdf) == 'rowname'))]
        write.table(mdf, paste0(header, '_', ident, '_', i, '_sig_snps_basic_data.tsv'), sep="\t")
    }
}

check.overlap.sig.snps <- function(idents, pvalue="emp_p_adj", threshold=0.05, fixed=FALSE, remx=FALSE) {
    if (fixed) {
        tail <- "fixed_"
    } else {
        tail <- ""
    }
    if (remx) {
        ltail <- "_remx"
    } else {
        ltail <- "_remx"
    }
    ident = (idents[1]-1)/3+1
    snp_info <- subset(all_df,  all_df[,'ident'] == ident & all_df[,"time"] == 1)
    sgene <- NULL
    sig_emp_adf <- NULL
    for (i in idents) {
        df <- read.table(paste0('fisher_stat_obs_emp_', tail, i, ltail, '.tsv'), header=T, sep="\t")
        if (is.null(sig_emp_adf)) sig_emp_adf <- cbind(df, time=i)
        else sig_emp_adf <- rbind(sig_emp_adf, cbind(df, time=i))
        sig_emp_df <- df[df[,pvalue] < threshold,]
        print(head(sig_emp_df))
        print(dim(sig_emp_df))
        if (is.null(sgene)) sgene <- cbind(sig_emp_df[,'gene'], rep((i-1)%%3, dim(sig_emp_df)[1]))
        else sgene <- rbind(sgene, cbind(sig_emp_df[,'gene'], rep((i-1)%%3, dim(sig_emp_df)[1])))
    }
    colnames(sgene) <- c('gene', 'time')
    print(sgene)
    occurrence <- sapply(unique(sgene[,'gene']), function(x) {
        return(length(which(x == sgene[,'gene'])))
    })
    print(length(occurrence[occurrence == 3]))
    print(length(occurrence[occurrence == 2]))
    print(length(occurrence[occurrence == 1]))
    gene_list = lapply(1:3, function(x){return(unique(sgene[,'gene'])[which(occurrence == x)])})
    saveRDS(gene_list, paste0('overlapped_sig_snps_', idents[1], tail, ltail, '.rds'))
    temp <- rdf[[ident]][[1]]
    for (x in gene_list[[3]]) {
        print(temp[temp[,'gene'] %in% x,])
    }
    plot.cov.and.ase(gene_list, (idents[1]-1)/3+1, paste0('overlapped_sig_snps_', idents[1], '_', tail, ltail)) # single model
    sig.snps.across.time(gene_list, sgene, sig_emp_adf, (idents[1]-1)/3+1, paste0('overlapped_sig_snps_', idents[1], '_', tail, ltail), ident, threshold=threshold) # paired model
}

plot.box.sig.only <- function(data, header, cor_label) {
    require(ggpubr)
    require(ggsci)
    require(viridis)
    for (norm in c('', '_weight')) {
        df  <- data.frame(data)
        if (norm != '') {
            df <- data %>% group_by_at(c('ident', 'ensemblID', 'significance', 'sig_1', 'sig_2', 'sig_3')) %>% summarise(
                    gene=min(gene),
                    first_pair=mean(first_pair, na.rm=TRUE),
                    second_pair=mean(second_pair, na.rm=TRUE)
            )
            df <- data.frame(df)
        }
        write.table(data, file=paste0(header, norm, '_', cor_label, '_data.tsv'), row.names=FALSE, quote=F, sep="\t")
        next
        print(head(df))
        print(sum(df[,'sig_1']))
        print(sum(df[,'sig_2']))
        print(sum(df[,'sig_3']))
        print(table(df[,'significance']))
        print(head(df[df[,'significance'] != 'others',]))
        all_data = NULL
        for (col in 1:2) {
            label = c('first', 'second')[col]
            print(colnames(df))
            sig = (df[,'sig_1'] & df[,paste0('sig_', col+1)])
            print(which(sig))
            temp = data.frame(time=label, significance='Not', correlation=df[!sig, paste0(label, '_pair')])
            if (length(which(sig)) > 0) {
                temp = rbind(temp, data.frame(time=label, significance='Deviated', correlation=df[sig, paste0(label, '_pair')]))
                print(wilcox.test(x=temp[temp[,'significance'] == 'Deviated', 'correlation'], y=temp[temp[,'significance'] == 'Not','correlation'], alternative='greater', paired=FALSE))
            }
            if (is.null(all_data)) all_data = temp
            else all_data = rbind(all_data, temp)
        }
        print('all_data')
        print(head(all_data))
        print(wilcox.test(x=all_data[all_data[,'significance'] == 'Deviated', 'correlation'], y=all_data[all_data[,'significance'] == 'Not','correlation'], alternative='greater', paired=FALSE))
        stopifnot(dim(all_data)[1] > 1)
        # g <- ggboxplot(all_data, x='time', y='correlation', fill='significance')+scale_fill_npg()
        g <- ggboxplot(all_data, x='time', y='correlation', fill='significance')+scale_fill_manual(values=c('white', 'gray'))
        # g <- facet(g, facet.by="time", scales='free_y')
        pdf(paste0('boxplot_', header, '_', cor_label, '_sig_only_all', norm, '.pdf'))
        plot(ggpar(g, font.main=c(20, "bold", "black"), font.x=16, font.y=16))
        dev.off()
        # g <- ggboxplot(all_data, x='significance', y='correlation', fill='significance')+scale_fill_npg()
        g <- ggboxplot(all_data, x='significance', y='correlation', fill='significance')+scale_fill_manual(values=c('white', 'gray'))
        # g <- facet(g, facet.by="time", scales='free_y')
        pdf(paste0('boxplot_', header, '_', cor_label, '_sig_only_all', norm, '_integrated.pdf'), width=4, height=6)
        plot(ggpar(g, font.main=c(20, "bold", "black"), font.x=16, font.y=16))
        dev.off()
    }
}


draw.sig.snps.cor <- function(fixed=TRUE, remx=TRUE) {
    if (fixed) {
        tail <- "fixed_"
    } else {
        tail <- ""
    }
    if (remx) {
        ltail <- '_remx'
    } else {
        ltail <- ''
    }
    for (cor_label in c('scor', 'cor')) {
        all_df <- NULL
        for (ident in seq(1, 13, 3)) {
            print(ident)
            header = paste0('overlapped_sig_snps_', ident, '_', tail, ltail)
            print(paste0(header, '_', cor_label, '_comp.csv'))
            df <- read.table(paste0(header, '_', cor_label, '_comp.csv'), sep=',', header=T)
            df <- cbind(df, ident=ident)
            if (is.null(all_df)) all_df <- df
            else all_df <- rbind(all_df, df)
        }
        print(head(all_df))
        plot.box.sig.only(all_df, paste0('overlapped_sig_snps_all__', tail, ltail), cor_label)
    }
}

if (!file.exists('black_list_home_snps.rds')) {
    extract.homozygous.snp()
}
fixed=FALSE
remx=TRUE
check.overlap.sig.snps(c(1:3), fixed=fixed, remx=remx)
check.overlap.sig.snps(c(4:6), fixed=fixed, remx=remx)
check.overlap.sig.snps(c(7:9), fixed=fixed, remx=remx)
check.overlap.sig.snps(c(10:12), fixed=fixed, remx=remx)
check.overlap.sig.snps(c(13:15), fixed=fixed, remx=remx)
draw.sig.snps.cor(fixed, remx=remx)
