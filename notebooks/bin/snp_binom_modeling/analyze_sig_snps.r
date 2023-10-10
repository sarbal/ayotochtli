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

print(head(rdf[[1]][[1]]))
print(head(sdf[[1]][[1]]))
print(head(sig_df))
print(head(sig_t0))
print(head(all_df))

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
    }
}

plot.accum.plot <- function() {
    return()
}
plot.mov.ave.plot <- function() {
    return()
}

plot.hist.group.correlation <- function(gene_list, ident, source, header) {
    require(reshape2)
    for (time in c(2, 3)) {
        df = source[source[,'time'] == time,]
        for (cor_label in c('scor', 'cor')) {
            tdf <- df[,c('gene', cor_label)]
            print(dim(tdf))
            tdf <- cbind(tdf, sapply(tdf[,'gene'], function(x){
                for (n in 1:length(gene_list)){
                    if (any(x == gene_list[[n]])) return(n)
                }
                return('others')
            }))
            colnames(tdf) <- c('gene', 'cor', 'significance')
            print(head(tdf))
            # g <- ggdensity(tdf, x = "cor",
            # add = "mean", color = "significance", fill = "significance")+
            g <- ggecdf(tdf, x = "cor",
            color = "significance")+xlim(1, -1)+
            scale_color_npg()+scale_fill_npg()
            png(paste0(header, '_', time, '_', cor_label, '.png'))
            plot(g)
            dev.off()
        }
    }
}
plot.cor.class <- function(gene_list, ident, source, header) {
    require(reshape2)
    for (cor_label in c('scor', 'cor')) {
        print(head(source))
        tdf <- source[,c('gene', cor_label, 'time')]
        tdf <- tdf[!duplicated(tdf[,c('gene', 'time')]),]
        print(dim(tdf))
        tdf <- cbind(tdf, sapply(tdf[,'gene'], function(x){
            for (n in 1:length(gene_list)){
                if (any(x == gene_list[[n]])) return(n)
            }
            return('others')
        }))
        for (i in 1:3) {
            tdf <- cbind(tdf, unlist(sapply(tdf[,'gene'], function(x) {
                return(any(x == gene_list[[i]]))
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
        # print(head(mdf))
        # print(head(all_df))
        # print(head(tdf))
        stopifnot(!is.na(mdf[1,'ensemblID']))
        write.table(mdf, file=paste0(header, '_', cor_label, '_comp.csv'), sep=',')
        plot.bar.rep.ratio(tdf, header, cor_label)
        tdf <- subset(tdf, tdf[,'significance'] != 'others')
        print(head(tdf))
        print(dim(tdf))
        g <- ggscatter(tdf, x = 'first_pair', y = 'second_pair', color = "significance")+
         scale_color_npg()
        png(paste0('scatter_', header, '_', cor_label, '.png'))
        plot(g)
        dev.off()
    }
    # detach(plyr)
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

convert.pn.label <- function(vec) {
    vec <- sapply(vec, function(x) {
        if (x == 'pos') return('Positive')
        else if (x == 'neg') return('Negative')
        else if (x == 'inc') return('-')
        return(NA)
    })
    return(vec)
}
add.rep.cor.each <- function(fixed_time, another_time, neg=FALSE) {
    return(unlist(sapply(1:length(fixed_time), function(i) {
        if (is.na(fixed_time[i]) || is.na(another_time[i])) return(NA)
        if ((fixed_time[i] > 0 && !neg) || (fixed_time[i] < 0 && neg)) {
            if (another_time[i] > 0) return('Positive')
            else if (another_time[i] < 0) return('Negative')
            else return('-')
        }
        return(NA)
    })))
}

plot.density.accum <- function(df, time, header, cor_label, norm) {
    print('start')
    print(head(df))
    g <- ggdensity(df, x = time,
        color = "black", fill = "significance", bw=3)+scale_fill_npg()
    pdf(paste0('density_', header, '_', cor_label, '_', time, '_all', norm, '.pdf'), width=7)
    plot(g)
    dev.off()
    g <- ggviolin(df, y = time, x = "significance",
        color = "black", fill = "significance", bw=3, add = c("mean_sd"))+scale_fill_npg()
    pdf(paste0('violin_', header, '_', cor_label, '_', time, '_all', norm, '.pdf'), width=7)
    plot(g)
    dev.off()
    g <- ggplot(df, aes_string(x = time, y= "..density..",  fill = "significance"))+geom_histogram(, position='dodge')+scale_fill_npg()+theme_bw()
    pdf(paste0('hist_', header, '_', cor_label, '_', time, '_all', norm, '.pdf'), width=10)
    plot(g)
    dev.off()
    df[,'count'] = 1
    df <- df[order(df[,time], decreasing=TRUE),]
    print(head(df))
    require(plyr)
    eval(parse(text=paste0('df.t <- ddply(df, .(significance), transform, cy = cumsum(count))')))
    detach("package:plyr", unload=TRUE)
    df.t <- data.frame(df.t)
    dict <- list()
    for (uni in unique(df.t[,'significance'])) {
        dict[[uni]] = max(df.t[df.t[,'significance'] == uni, 'cy'])
    }
    print(dict)
    print('here')
    print('death')
    print(dict)
    # print(head(df.t))
    # print(head(df.t[, 'significance']))
    df.t[,'cy'] <- sapply(1:dim(df.t)[1], function(x) {return(df.t[x, 'cy']/dict[[df.t[x, 'significance']]])})
    print(head(df.t))
    df.t[,time] <- as.numeric(df.t[,time])
    print('?????')
    # print(colnames(df.t))
    # print(df.t[,'significance'])
    df.t <- data.frame(df.t)
    g <- ggplot(df.t, aes_string(x=time, y="cy"))+geom_line(aes(color=significance), size=2, alpha=0.9)+scale_color_npg()+theme_bw()+xlim(1, -1)
    pdf(paste0('line_', header, '_', cor_label, '_', time, '_all', norm, '.pdf'))
    # png(paste0('line_', header, '_', cor_label, '_', time, '_all', norm, '.png'))
    plot(g)
    dev.off()
    pdf(paste0('line_', header, '_', cor_label, '_', time, '_all', norm, '.pdf'))
    plot(g)
    dev.off()
}

plot.density.sig <- function(df, time, header, cor_label, norm) {
    
    print('start')
    print(head(df))
    g <- ggdensity(df, x = time,
        color = "black", fill = "significance", bw=3)+scale_fill_npg()
    pdf(paste0('density_', header, '_', cor_label, '_', time, '_all', norm, '.pdf'), width=7)
    plot(g)
    dev.off()
    g <- ggviolin(df, y = time, x = "significance",
        color = "black", fill = "significance", bw=3, add = c("mean_sd"))+scale_fill_npg()
    pdf(paste0('violin_', header, '_', cor_label, '_', time, '_all', norm, '.pdf'), width=7)
    plot(g)
    dev.off()
    g <- ggplot(df, aes_string(x = time, y= "..density..",  fill = "significance"))+geom_histogram(, position='dodge')+scale_fill_npg()+theme_bw()
    pdf(paste0('hist_', header, '_', cor_label, '_', time, '_all', norm, '.pdf'), width=10)
    plot(g)
    dev.off()
    df[,'count'] = 1
    df <- df[order(df[,time], decreasing=TRUE),]
    print(head(df))
    require(plyr)
    eval(parse(text=paste0('df.t <- ddply(df, .(significance), transform, cy = cumsum(count))')))
    detach("package:plyr", unload=TRUE)
    df.t <- data.frame(df.t)
    dict <- list()
    for (uni in unique(df.t[,'significance'])) {
        dict[[uni]] = max(df.t[df.t[,'significance'] == uni, 'cy'])
    }
    print(dict)
    print('here')
    print('death')
    print(dict)
    # print(head(df.t))
    # print(head(df.t[, 'significance']))
    df.t[,'cy'] <- sapply(1:dim(df.t)[1], function(x) {return(df.t[x, 'cy']/dict[[df.t[x, 'significance']]])})
    print(head(df.t))
    df.t[,time] <- as.numeric(df.t[,time])
    print('?????')
    # print(colnames(df.t))
    # print(df.t[,'significance'])
    df.t <- data.frame(df.t)
    g <- ggplot(df.t, aes_string(x=time, y="cy"))+geom_line(aes(color=significance), size=2, alpha=0.9)+scale_color_npg()+theme_bw()+xlim(1, -1)
    pdf(paste0('line_', header, '_', cor_label, '_', time, '_all', norm, '.pdf'))
    # png(paste0('line_', header, '_', cor_label, '_', time, '_all', norm, '.png'))
    plot(g)
    dev.off()
    pdf(paste0('line_', header, '_', cor_label, '_', time, '_all', norm, '.pdf'))
    plot(g)
    dev.off()
}
plot.box.sig.only <- function(df, header, cor_label) {
    require(ggpubr)
    require(ggsci)
    require(viridis)
    for (norm in c('', '_weight')) {
        temp  <- data.frame(df)
        if (norm != '') {
            temp <- df %>% group_by_at(c('ident', 'ensemblID', 'significance', 'sig_1', 'sig_2', 'sig_3')) %>% summarise(
                    gene=min(gene),
                    first_pair=mean(first_pair, na.rm=TRUE),
                    second_pair=mean(second_pair, na.rm=TRUE)
            )
            temp <- data.frame(temp)
        }
        print(colnames(temp))
        print(head(temp))
        temp[,'significance'] <- factor(temp[,'significance'], levels=c('Significant (n=3)', 'Significant (n=2)', 'Significant (n=1)', 'Non-significant'))
        print(norm)
        print(head(temp))
        plot.density.sig(temp, time, header, cor_label, norm)        
    }
}


plot.box.rep.ratio <- function(df, header, cor_label) {
    require(ggpubr)
    require(ggsci)
    require(viridis)
    df[,'significance'] <- unlist(sapply(df[,'significance'], function(x) {
        if (x == '1') return('Significant (n=1)')
        if (x == '2') return('Significant (n=2)')
        if (x == '3') return('Significant (n=3)')
        if (x == 'others') return('Non-significant')
    }))
    # df <- subset(df[df[,'ident'] == 1,])
    for (norm in c('', '_weight')) {
        temp  <- data.frame(df)
        if (norm != '') {
            temp <- df %>% group_by_at(c('ident', 'ensemblID', 'significance')) %>% summarise(
                    gene=min(gene),
                    first_pair=mean(first_pair, na.rm=TRUE),
                    second_pair=mean(second_pair, na.rm=TRUE)
            )
            temp <- data.frame(temp)
        }
        print(colnames(temp))
        print(head(temp))
        temp[,'significance'] <- factor(temp[,'significance'], levels=c('Significant (n=3)', 'Significant (n=2)', 'Significant (n=1)', 'Non-significant'))
        print(norm)
        print(head(temp))
        if (TRUE) {
            for (time in c('first_pair', 'second_pair')) {
                plot.density.accum(temp, time, header, cor_label, norm)
            }
        } else {
            print(c('norm', norm))
            print(head(temp))
            print('set labels')
            print('rep')
            temp[,'rep'] <- add.rep.cor.label(temp, TRUE)
            print('first')
            temp[,'first_pos'] <- add.rep.cor.each(temp[,'first_pair'], temp[,'second_pair'])
            temp[,'first_neg'] <- add.rep.cor.each(temp[,'first_pair'], temp[,'second_pair'], neg=TRUE)
            print('second')
            temp[,'second_pos'] <- add.rep.cor.each(temp[,'second_pair'], temp[,'first_pair'])
            temp[,'second_neg'] <- add.rep.cor.each(temp[,'second_pair'], temp[,'first_pair'], neg=TRUE)
            for (col in c('rep', 'first_pos', 'second_pos', 'first_neg', 'second_neg')) {
                print(col)
                # print(head(ttemp))
                if (col == 'rep') {
                    temp[,col] <- factor(temp[,col], level=c('Positive', 'Negative', '-'))
                } else {
                    temp[,col] <- factor(temp[,col], level=c('Positive', 'Negative'))
                }
                ttemp <- temp[temp[,col] != 'nan',]
                ttemp <- ttemp[!is.na(ttemp[,col]),]
                print(head(ttemp))
                sum_tdf <- ttemp %>% group_by_at(c('ident', col, 'significance')) %>% summarise(
                    count=n()
                )
                sum_tdf <- data.frame(sum_tdf)
                g <- ggboxplot(sum_tdf, x=col, y='count', fill=col, add = "jitter", outlier.shape = NA)+scale_fill_npg()
                g <- facet(g, facet.by="significance", scales='free_y')
                pdf(paste0('boxplot_', header, '_', cor_label, '_', col, '_all', norm, '.pdf'))
                plot(ggpar(g, font.main=c(20, "bold", "black"), font.x=16, font.y=16))
                dev.off()
            }
        }
    }
}

plot.bar.rep.ratio <- function(df, header, cor_label) {
    tdf <- df[!is.na(df[,'rep']),]
    tdf[,'rep'] <- factor(tdf[,'rep'], level=c('pos', 'neg', 'inc', 'nan'))
    adf <- subset(tdf, tdf[,'significance'] != 'others')
    stat = 'bin'
    if (cor_label == 'scor') stat = 'count'
    g <- ggplot(adf, aes(significance, fill=rep) )+geom_bar(position='stack', stat='count')+scale_fill_npg()
    png(paste0('barplot_', header, '_', cor_label, '_sig.png'))
    plot(g)
    dev.off()
    adf <- subset(tdf, tdf[,'significance'] == 'others')
    g <- ggplot(adf, aes(x=significance, y = (..count..)/sum(..count..), fill=rep) )+geom_bar(position='stack')+scale_fill_npg()
    png(paste0('barplot_', header, '_', cor_label, '_all.png'))
    plot(g)
    dev.off()
}

plot.cor.class.all <- function(fixed=TRUE, remx=TRUE) {
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
        # plot.bar.rep.ratio(all_df, paste0('overlapped_sig_snps_all__', tail, ltail), cor_label)
        # plot.box.rep.ratio(all_df, paste0('overlapped_sig_snps_all__', tail, ltail), cor_label)
        plot.box.sig.only(all_df, paste0('overlapped_sig_snps_all__', tail, ltail), cor_label)
    }
}

plot.sig.snps.cor <- function(gene_list, emp_p_df, pvalue, header, ident, threshold=0.05) {
    require(ggplot2)
    require(ggsci)
    require(ggpubr)
    print(ident)

    print(head(sdf[[ident]][[1]]))
    source = cbind(sdf[[ident]][[1]], time=2)
    source = rbind(source, cbind(sdf[[ident]][[2]], time=3))
    source[,'time'] = factor(source[,'time'])
    # plot.accum.plot(pvalue)
    # plot.mov.ave.plot()
    # plot.hist.group.correlation(gene_list, ident, source, header)
    plot.cor.class(gene_list, ident, source, header)
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
    # plot.cov.and.ase(gene_list, (idents[1]-1)/3+1, paste0('overlapped_sig_snps_', idents[1], '_', tail, ltail)) # single model
    plot.sig.snps.cor(gene_list, sig_emp_adf, (idents[1]-1)/3+1, paste0('overlapped_sig_snps_', idents[1], '_', tail, ltail), ident, threshold=threshold) # paired model
    # sdf[[ident]][[1]]
}

# extract.homozygous.snp()
fixed=FALSE
remx=TRUE
# check.overlap.sig.snps(c(1:3), fixed=fixed, remx=remx)
# check.overlap.sig.snps(c(4:6), fixed=fixed, remx=remx)
# check.overlap.sig.snps(c(7:9), fixed=fixed, remx=remx)
# check.overlap.sig.snps(c(10:12), fixed=fixed, remx=remx)
# check.overlap.sig.snps(c(13:15), fixed=fixed, remx=remx)
plot.cor.class.all(fixed, remx=remx)


