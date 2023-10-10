require(ggplot2)
require(ggpubr)
require(reshape2)
require(viridis)
require(parallel)
require(dplyr)
sdf <- readRDS('simple_independent_result.rds')
# sdf <- readRDS('simple_t0_result.rds')
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
print(cov)
print(median(cov))
all_data <- log10(all_data)
q <- quantile(all_data, c(25, 50, 75, 90, 95)/100)
all_data <- data.frame(cov=all_data)
print(q)
g <- ggdensity(all_data, x='cov', rug = TRUE)+geom_vline(xintercept=q[1], linetype='dotted')+geom_vline(xintercept=q[2], linetype='dotted')+geom_vline(xintercept=q[3], linetype='dotted')

pdf(paste0('density_coverage_all_quads.pdf'))
plot(g)
dev.off()

