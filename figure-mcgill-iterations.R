library(ggplot2)
library(data.table)
iterations.dt <- data.table::fread("figure-mcgill-iterations-data.csv")
bound.dt <- fread("figure-mcgill-iterations-bound.csv")
iterations.dt[order(N.data)]
range(iterations.dt$N.data)

ggplot()+
  geom_point(aes(
    N.data, sum.splits, color=max.segments),
    data=iterations.dt)+
  scale_x_log10()+
  scale_y_log10()

size <- 1
gg <- ggplot()+
  geom_point(aes(
    N.data, sum.splits),
    shape=".",
    color="grey50",
    data=iterations.dt)+
  geom_line(aes(
    N.data, sum.splits, color=case, size=case),
    data=bound.dt)+
  facet_grid(. ~ max.segments, labeller=label_both)+
  scale_x_log10(
    "Number of data to segment",
    labels=scales::label_log())+
  scale_y_log10(
    "Number of candidate\nsplits to consider",
    labels=scales::label_log())+
  theme_bw()+
  theme(
    axis.text=element_text(size=12),
    panel.spacing=grid::unit(0,"lines"))+
  coord_cartesian(ylim=c(1e2, 1e7))+
  scale_size_manual(values=c(
    best=1.5, best.heuristic=0.75, worst=0.75))
png("figure-mcgill-iterations.png", width=7, height=2, units="in", res=200)
print(gg)
dev.off()

iterations.dt[, round.N.data := 10^(round(log10(N.data)))]
iterations.dt[, num.sum.splits := as.numeric(sum.splits)]
iterations.stats <- dcast(
  iterations.dt,
  round.N.data + max.segments ~ .,
  fun.aggregate=list(median, min, max, length),
  value.var="num.sum.splits")
iterations.stats[, case := "empirical"]
case.colors <- c(
  empirical="black",
  best="#E41A1C",
  worst="#377EB8",
  best.heuristic="#4DAF4A")#, "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF")
ggplot()+
  scale_color_manual(values=case.colors)+
  scale_fill_manual(values=case.colors, guide="none")+
  geom_ribbon(aes(
    round.N.data,
    ymin=num.sum.splits_min,
    ymax=num.sum.splits_max,
    fill=case),
    alpha=0.5,
    data=iterations.stats)+
  geom_line(aes(
    round.N.data, num.sum.splits_median, color=case),
    size=size,
    data=iterations.stats)+
  facet_grid(. ~ max.segments, labeller=label_both)+
  scale_x_log10(
    "Number of data to segment")+
  scale_y_log10(
    "Number of candidate splits to consider")+
  geom_line(aes(
    N.data, sum.splits, color=case),
    size=size,
    data=bound.dt)+
  theme_bw()+
  theme(panel.spacing=grid::unit(0,"lines"))+
  coord_cartesian(ylim=c(1e2, 1e7))
