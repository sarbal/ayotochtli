
  1
  2
  3
  4
  5
  6
  7
  8
  9
 10
 11
 12
 13
 14
 15
 16
 17
 18
 19
 20
 21
 22
 23
 24
 25
 26
 27
 28
 29
 30
 31
 32
 33
 34
 35
 36
 37
 38
 39
 40
 41
 42
 43
 44
 45
 46
 47
 48
 49
 50
 51
 52
 53
 54
 55
 56
 57
 58
 59
 60
 61
 62
 63
 64
 65
 66
 67
 68
 69
 70
 71
 72
 73
 74
 75
 76
 77
 78
 79
 80
 81
 82
 83
 84
 85
 86
 87
 88
 89
 90
 91
 92
 93
 94
 95
 96
 97
 98
 99
100
101
102
103
104
105
106
107
108
109
110
111
112
113
114
115
116
117
118
119
120
121
122
123
124
125
126
127
128
129
130
131
132
133
134
135
136
137
138
139
140
141
142
143
144
145
146
147
148
149
150
151
152
153
154
155
156
157
158
159
160
161
162
163
164
165
166
167
168
169
170
171
172
173
174
175
176
177
178
179
180
181
182
183
184
185
186
187
188
189
190
191
192
require(ggplot2)
require(ggpubr)
require(reshape2)
require(viridis)
require(parallel)
sdf <- readRDS('rds_files/simple_t0_result.rds')

read.fisher.data <- function(fix, i, remx, label) {
    fname = paste0("table/fisher_stat_obs_emp", fix, i, remx, '.tsv')
    if (!file.exists(fname)) {
        return(NULL)
    }
    a <- read.table(fname, header=T, sep="\t")
    ident = floor((i-1)/3+1)
    time = ((i-1)%%3+1)
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

merge.fisher.correlation <- function(cor, fix, label, remx) {
    pdata = NULL
    for (i in 1:15) { 
        result <- read.fisher.data(fix, i, remx, label)
        if (is.null(result)) next
        # size <- c(size, result[[1]])
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
        if (is.null(pdata)) pdata <- rtemp
        else pdata <- rbind(pdata, rtemp)
    }
    pdata <- cbind(pdata, snp=rep(remx, dim(pdata)[1]))
    return(pdata)
}

plot.hist.all.snp.groups <- function() {
    for (cor in c('cor')) {
        for (fix in c('_')) {
            for (label in c('emp_p_adj')) {
                all_table <- NULL
                for (remx in c('', '_remx', '_extx')) {
                    pdata <- merge.fisher.correlation(cor, fix, label, remx)
                    if (is.null(all_table)) all_table <- pdata
                    else all_table <- rbind(all_table, pdata)
                }
                all_table <- all_table[!is.na(all_table[,'cor']),] # remove first time point                           
                all_table[,4] <- sapply(all_table[,4], function(x){if(x == "") return("All"); if(x == "_remx") return("No X"); return("X");})
                df <- all_table %>%   
                        group_by(snp) %>% 
                        mutate(rank = rank(value, ties.method = "first"))
                df <- data.frame(df)
                colnames(df)[4] = "SNP"
                plot.rank.cor.binned(cor, fix, label, df)
                plot.cor.significance(cor, fix, label, df)
            }
        }
    }
}


compute.sum.count <- function(tdf) {
    tdf <- cbind(tdf, significance=unlist(sapply(tdf[,'value'], function(x){if(x >= -log10(0.05)) return('Significant'); return('Non significant')})))
    tdf <- tdf[order(tdf[,'cor'], decreasing=TRUE),]
    tdf <- cbind(tdf, count=rep(1, dim(tdf)[1]))
    mdf <- tdf %>% group_by(crank, significance) %>%
        mutate(sum_count = cumsum(count))
    mdf <- mdf[!duplicated(mdf[,c('crank', 'SNP', 'significance')], fromLast = TRUE),]
    mdf <- data.frame(mdf)
    mdf <- mdf[order(mdf[,'crank'], decreasing=TRUE),]
    mdf <- mdf %>% group_by(significance) %>% mutate(cum_count = cumsum(sum_count))
    return(mdf)
}

compute.cum.count <- function(mdf) {
    wide = mdf[,c('crank', 'SNP', 'significance', 'cum_count')] %>% spread(significance, cum_count)
    wide = data.frame(wide)
    for (i in 1:(dim(wide)[1]-1)) {
        for (j in 3:4) {
            if (is.na(wide[i,j])) {
                wide[i,j] = max(wide[(i+1):dim(wide)[1],j], na.rm=TRUE)
            }
        }
    }
    wide = cbind(wide, ratio=(wide[,4]+1)/(wide[,4]+wide[,3]+1))
    wide[,5] = wide[,5]/wide[1,5]
    wide[,1] <- as.numeric(as.character(wide[,1]))
    return(wide)
}

plot.cor.significance <- function(cor, fix, label, df) {
    require(dplyr)
    require(tidyverse)
    for (mode in c('All', 'No X', 'X')) {
        colors <- list()
        if (mode == 'No X') {
            colors[['Significant']] <- rgb(91/255, 89/255, 166/255)
            tdf <- df[df[,'SNP'] == 'No X',]
        } else if (mode == 'All') {
            colors[['Significant']] <- rgb(20/255, 175/255, 158/255)
            tdf <- df[df[,'SNP'] == 'All',]
        } else {
            colors[['Significant']] <- rgb(184/255, 82/255, 157/255)
            tdf <- df[df[,'SNP'] == 'X',]
        }
        colors[['Non significant']] <- "gray"
        if (cor == 'scor') {
            xlab = "Spearman Correlation Coefficient"
        } else {
            xlab = "Pearson Correlation Coefficient"
        }
        ylab = "Cumulative ratio of deviated SNPs"
        if (cor == 'cor') {
            tdf <- cbind(tdf, crank=cut(tdf[,'cor'], seq(-1.0, 1, 0.1), labels=seq(-1.0, 0.9, 0.1), include.lowest=TRUE))
        } else {
            tdf <- cbind(tdf, crank=cut(tdf[,'cor'], seq(-1.0, 1, 0.1), labels=seq(-1.0, 0.9, 0.1), include.lowest=TRUE))
        }
        mdf <- compute.sum.count(tdf)
        write.table(mdf, paste0('sum_', fix, '_', label, '_', cor, '_', mode, '.tsv'))
        wide <- compute.cum.count(mdf)
        write.table(wide, paste0('cumsum_', fix, '_', label, '_', cor, '_', mode, '.tsv'))
        p <- ggplot(wide, aes(x=crank, y=ratio)) + geom_point() + geom_smooth(method = "loess", se=FALSE, color=colors[[1]]) +
                theme_bw()+theme(text = element_text(size=20)) +  xlab(xlab)+ylab(ylab) 
        pdf(paste0('plot_cor_ratio_', fix, '_all', '_', label, '_', cor, '_', mode, '.pdf'))
        plot(p)
        dev.off()
    }
}


plot.rank.cor.binned <- function(cor, fix, label, df) {
    count <- table(df[,'SNP'])
    df[,'rank'] <- unlist(mclapply(1:dim(df)[1], function(x){return(df[x,'rank']/count[names(count) == df[x,'SNP']])}, mc.cores=10))
    write.table(df, paste0('rank_data_', fix, '_all_', label, '_', cor, '.csv'), sep=",")
    methods <- c('lm')
    for (mode in c('All', 'No X')) {
        colors <- list()
        if (mode == 'No X') {
            colors[['No X']] <- rgb(91/255, 89/255, 166/255)
            tdf <- df[df[,'SNP'] == 'No X',]
        } else if (mode == 'All') {
            colors[['All']] <- rgb(20/255, 175/255, 158/255)
            tdf <- df[df[,'SNP'] == 'All',]
        }
        tdf <- rbind(tdf, df[df[,'SNP'] == 'X',])
        colors[['X']] <- rgb(184/255, 82/255, 157/255)
        ylab = "Pearson Correlation Coefficient"
        if (cor == 'scor') ylab = "Spearman Correlation Coefficient"
        xlab = "Ranking order of -log10 FDR"
        tdf[,'qrank'] <- cut(tdf[,'rank'], seq(0, 1, 0.2), labels=seq(0.2, 1.0, 0.2), include.lowest=TRUE)
        label = c('0.2'='0.2', '0.4'='0.4', '0.6'='0.6', '0.8'='0.8', '1'='1.0')
        tdf[,'fisher'] <- tdf[,'value']
        tdf[,'value'] <- cut(tdf[,'value'], seq(0, 3, 1), labels=seq(1, 3, 1), include.lowest=TRUE)
        label = c('1'='0.1', '2'='0.01', '3'='0.001')
        if (mode == 'All')
            write.table(tdf, file='test_output.txt')
        p <- ggplot(subset(tdf[tdf[,'SNP'] == 'X',])) + geom_boxplot(aes(value, cor), fill=colors[['X']])+theme_bw()+theme(text=element_text(size=20))+xlab(xlab)+ylab(ylab)
        p <- p+scale_x_discrete(labels=label)
        pdf(paste0('fdrplot_binned_box_p_', fix, '_all', '_', label, '_', cor, '_X.pdf'))
        plot(p)
        dev.off()
        p <- ggplot(subset(tdf[tdf[,'SNP'] == mode,])) + geom_boxplot(aes(value, cor), fill=colors[[mode]])+theme_bw()+theme(text=element_text(size=20))+xlab(xlab)+ylab(ylab)
        p <- p+scale_x_discrete(labels=label)
        pdf(paste0('fdrplot_binned_box_p_', fix, '_all', '_', label, '_', cor, '_', mode, '.pdf'))
        plot(p)
        dev.off()
    }
}

plot.hist.all.snp.groups()
