source("packages.R")
N.data.exp <- 6
N.data <- 2^N.data.exp
print(N.data)
set.seed(1)
data.list <- list(
  rpois1000=rpois(N.data,1000),
  linear=1:N.data,
  zero.ten=rep(c(0,1,10,11),l=N.data),
  zero.one=rep(c(0,1),l=N.data))
out.dt.list <- list()
data.name <- "zero.one"
loss <- "meanvar_norm"
for(loss in c("l1", "poisson", "meanvar_norm", "mean_norm")){
  for(data.name in names(data.list)){
    data.vec <- data.list[[data.name]]
    binseg.fit <- binsegRcpp::binseg(loss, data.vec)
    clist <- binsegRcpp::get_complexity(binseg.fit)
    clist$iterations[case=="worst"]
    for(out.name in names(clist)){
      out.dt.list[[out.name]][[paste(loss, data.name)]] <- data.table(
        loss, data.name, clist[[out.name]])
    }
  }
}
out <- lapply(out.dt.list, function(L)do.call(rbind, L))
out.wide <- dcast(
  out$iterations, 
  loss + data.name + segments ~ case, 
  value.var="cum.splits")
out.wide[empirical < best]

gg <- ggplot()+
  facet_grid(loss ~ data.name, labeller=label_both)+
  geom_line(aes(
    segments, cum.splits, color=case, size=case),
    data=out$iterations[case!="empirical"])+
  geom_point(aes(
    segments, cum.splits, color=case),
    shape=1,
    data=out$iterations[case=="empirical"])+
  scale_color_manual(
    values=binsegRcpp::case.colors,
    breaks=names(binsegRcpp::case.colors))+
  scale_x_log10()+
  scale_y_log10()+
  scale_size_manual(
    values=binsegRcpp::case.sizes,
    guide="none")+
  theme_bw()+
  theme(panel.spacing=grid::unit(0,"lines"))
print(gg)
png("figure-splits-loss-cum.png", width=14, height=6, units="in", res=200)
print(gg)
dev.off()

gg <- ggplot()+
  facet_grid(loss ~ data.name, labeller=label_both)+
  geom_line(aes(
    segments, splits, color=case, size=case),
    data=out$iterations[case!="empirical"])+
  geom_point(aes(
    segments, splits, color=case),
    shape=1,
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
print(gg)
png("figure-splits-loss.png", width=14, height=6, units="in", res=200)
print(gg)
dev.off()
