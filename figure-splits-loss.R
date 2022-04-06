source("packages.R")

N.data.exp <- 4
N.data <- 2^N.data.exp
print(N.data)
set.seed(1)
data.list <- list(
  flat=rnorm(N.data),
  linear=1:N.data,
  zero.one=rep(c(0,1),l=N.data))
out.dt.list <- list()
for(loss in c("l1", "poisson", "meanvar_norm", "mean_norm")){
  for(data.name in names(data.list)){
    data.vec <- data.list[[data.name]]
    binseg.fit <- binsegRcpp::binseg(loss, data.vec)
    clist <- binsegRcpp::get_complexity(binseg.fit)
    for(out.name in names(clist)){
      out.dt.list[[out.name]][[paste(loss, data.name)]] <- data.table(
        loss, data.name, clist[[out.name]])
    }
  }
}
out <- lapply(out.dt.list, function(L)do.call(rbind, L))
library(ggplot2)
gg <- ggplot()+
  facet_grid(loss ~ data.name, labeller=label_both)+
  geom_line(aes(
    segments, splits, color=case, size=case),
    data=out$iterations[case!="empirical"])+
  geom_point(aes(
    segments, splits, color=case),
    data=out$iterations[case=="empirical"])+
  geom_text(aes(
    x, y,
    label=label,
    color=case),
    hjust=1,
    data=out$totals)+
  scale_color_manual(
    values=binsegRcpp::case.colors,
    guide="none")+
  scale_size_manual(
    values=binsegRcpp::case.sizes,
    guide="none")+
  theme_bw()+
  theme(panel.spacing=grid::unit(0,"lines"))
png("figure-splits-loss.png", width=10, height=6, units="in", res=200)
print(gg)
dev.off()
