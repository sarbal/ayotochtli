require(ComplexHeatmap)
require(ggplot2)
require(reshape2)
require(tidyr)
require(viridis)
require(ggsci)
require(dplyr)
cov_quantile <- c(1.051153, 1.295567, 1.625827, 1.977724, 2.202761)
cov_vec <- ceiling(10**(cov_quantile))
print(cov_vec)

read.cov.dist <- function() {
    sdf <- readRDS('simple_independent_result.rds')
    cov <- c()
    all_data <- c()
    for (i in 1:length(sdf)) {
        for (j in 1:length(sdf[[i]])) {
            a <- sdf[[i]][[j]] %>% group_by(gene) %>% summarise(
                mean_cov = mean(cov)
            )
            a <- data.frame(a)
            print(head(a))
            all_data <- c(all_data, a[,'mean_cov'])
            cov <- c(cov, median(a[,'mean_cov']))
        }
    }
    return(all_data)
}
plot.quantile.value <- function(part, prob) {
    header = paste0('plot_quantile_cov_', prob, '.pdf')
    # print(part)
    part <- melt(as.matrix(part))
    colnames(part) <- c('cell', 'cov', 'AUROC')
    print(head(part))
    part[,'cell'] <- factor(part[,'cell'])
    part[,'cov'] <- factor(part[,'cov'])
    g <- ggplot(part, aes(x=cell, y=AUROC, group=cov, color=cov))+geom_line()+geom_point()+scale_color_npg()+theme_bw()
    pdf(header)
    plot(g)
    dev.off()
}

plot.real.cov.dist <- function(m, prob) {
    ind = 6
    print(head(m))
    all_data <- read.cov.dist()
    df <- data.frame()
    for (ind in c(5, 6, 7)) {
        auroc <- unlist(sapply(all_data, function(x){return(m[ind,which(colnames(m) == x)])}))
        df <- rbind(df, data.frame(AUROC=auroc, cell=rep(rownames(m)[ind], length(auroc))))
    }
    df[,'cell'] <- factor(df[,'cell'])
    pdf(paste0('plot_real_cov_dist_', prob, '.pdf'))
    g <- ggplot(df, aes(x=AUROC, group=cell, fill=cell))+geom_histogram(color='black', alpha=0.4, position='identity', bins=25)+xlim(0.5, 1)+theme_bw()+scale_color_npg()
    plot(g)
    dev.off()
}

for (prob in c(0.001, 0.01, 0.05)) {
    rep <- 100
    a <- read.table(paste0('output_auc_', rep, '_', prob, '_lm.txt'), header=T, sep=" ")
    colnames(a) <- c('cell', 'cov', 'quad', 'value')
    a <- subset(a, a[,'quad'] == 5)
    rownames(a) <- 1:dim(a)[1]
    a <- a[,c(1, 2, 4)]
    a[,'value'] <- as.numeric(a[,'value'])
    m <- dcast(a, cell ~ cov, value.var="value",  mean)
    rownames(m) <- m[,1]
    m <- m[,2:dim(m)[2]]
    hm <- Heatmap(m, cluster_rows=FALSE, cluster_columns=FALSE, col=viridis(50), show_column_names = FALSE)
    pdf(paste0('heatmap_auc_', prob, '.pdf'))
    draw(hm)
    dev.off()
    part <- m[,colnames(m) %in% cov_vec]
    plot.quantile.value(part, paste0(rep, '_', prob))
    plot.real.cov.dist(m, paste0(rep, '_', prob))


}
