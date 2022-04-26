library(data.table)
library(ggplot2)
iterations.dt <- data.table::fread("figure-neuroblastoma-iterations-data.csv")
bound.dt <- data.table::fread("figure-neuroblastoma-iterations-bounds.csv")
bound.wide <- dcast(
  bound.dt,
  N.data + max.segments ~ case,
  value.var = "sum.splits")
bound.and.data <- bound.wide[
  iterations.dt, on=.(N.data, max.segments), nomatch=0L]
## These three should be empty.
bound.and.data[sum.splits < best]
bound.and.data[best.heuristic < best]
bound.and.data[sum.splits > worst]
## These are real data sets for which the heuristic does not work.
(hilite <- bound.and.data[sum.splits < best.heuristic, .(
  N.data, max.segments, best, sum.splits, best.heuristic, worst)])
bound.and.data[sum.splits == best]
unique(
  bound.and.data[max.segments==10][, .(best.heuristic,best)]
)[best.heuristic == best]

ggplot()+
  geom_point(aes(
    N.data, sum.splits, color=max.segments),
    data=iterations.dt)+
  scale_x_log10()+
  scale_y_log10()

gg <- ggplot()+
  geom_line(aes(
    N.data, sum.splits, color=case, size=case),
    data=bound.dt)+
  geom_point(aes(
    N.data, sum.splits),
    shape=".",
    color="grey50",
    data=iterations.dt)+
  geom_point(aes(
    N.data, sum.splits),
    shape=1,
    data=hilite)+
  facet_grid(. ~ max.segments, labeller=label_both)+
  scale_x_log10(
    "Number of data to segment")+
  scale_y_log10(
    "Number of candidate\nsplits to consider")+
  theme_bw()+
  theme(panel.spacing=grid::unit(0,"lines"))+
  coord_cartesian(ylim=c(3e1, 2e5))+
  scale_size_manual(values=c(
    best=1.5, best.heuristic=0.75, worst=0.75))
png("figure-neuroblastoma-iterations.png", width=7, height=2, units="in", res=200)
print(gg)
dev.off()

