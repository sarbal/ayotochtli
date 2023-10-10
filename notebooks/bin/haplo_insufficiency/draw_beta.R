all_df = NULL

params_a = c(1, 2, 4, 10, 15, 30, 1, 1, 1, 4, 10, 30)
params_b = c(1, 2, 4, 10, 15, 30, 4, 10, 30, 1, 1, 1)
for (count in 1:length(params_a)) {
    x = seq(0, 1.01, 0.01)
    y = unlist(sapply(x, function(q){return(pbeta(q, params_a[count], params_b[count], ncp = 0, lower.tail = TRUE, log.p = FALSE))}))
    if (is.null(all_df)) {
        all_df = data.frame(x=x, y=y, a=params_a[count], b=params_b[count], count=count)
    } else {
        all_df = rbind(all_df, data.frame(x=x, y=y, a=params_a[count], b=params_b[count], count=count))   
    }

}
all_df[,'group'] = paste0(all_df[,'count'], '_', all_df[,'a'], '_', all_df[,'b'])
require(ggplot2)
require(ggsci)
require(viridis)
g <- ggplot(all_df, aes(x=x, y=y, group=group, color=group))+geom_line()+scale_color_viridis(discrete=TRUE)+theme_bw()
png("beta_sigmoid.png")
plot(g)
dev.off()
for (i in 10:12) {
    b <- all_df[all_df[,'count'] == i,]
    g <- ggplot(b, aes(x=x, y=y))+geom_line(lwd=3)+theme_bw()
    png(paste0("beta_sigmoid_", i, ".png"))
    plot(g)
    dev.off()

}