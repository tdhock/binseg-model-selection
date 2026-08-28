library(ggplot2)
library(data.table)
ares <- readRDS("figure-splits-time-data.rds")
ares$measurements <- ares$measurements[N>50]
ares$unit.col.vec <- c(seconds="median", "splits")
aref <- atime::references_best(ares)
aref$plot.references <- aref$references[fun.name %in% c("N log N", "N^2")]

plot(ares)

plot(aref)+scale_x_log10(
  limits=c(1, NA))+
  facet_grid(unit ~ expr.name, scales="free_y")

select.dt <- rowwiseDT(
  fun.name=, case=,
  "N log N", "best",
  "N^2", "worst")
aref$plot.references <- aref$references[select.dt, on=names(select.dt)]
plot(aref)+scale_x_log10(
  limits=c(1, NA))+
  facet_grid(unit ~ ., scales="free")

apred <- predict(aref, seconds=ares$seconds.limit, splits=1e6)
plot(apred)+
  geom_line(aes(
    N, reference, group=expr.name),
    data=aref$plot.references)

meas <- apred[["measurements"]][unit %in% apred$prediction$unit]
pred <- apred[["prediction"]]
one <- pred[, .SD[1], by=unit]
gg <- ggplot2::ggplot()+
  ggplot2::theme_bw()+
  ggplot2::facet_grid(unit ~ ., scales="free")+
  ggplot2::geom_hline(ggplot2::aes(
    yintercept=unit.value),
    color="grey",
    data=one)+
  ggplot2::geom_text(ggplot2::aes(
    100, unit.value, label=paste0(unit,"=",unit.value)),
    color="grey50",
    hjust=0,
    vjust=1.2, 
    data=one)+
  ggplot2::geom_ribbon(ggplot2::aes(
    N, ymin=min, ymax=max, fill=expr.name),
    data=meas[unit=="seconds"],
    alpha=0.5)+
  ggplot2::geom_line(ggplot2::aes(
    N, empirical, color=expr.name),
    size=1,
    data=meas)+
  geom_line(aes(
    N, reference, group=expr.name),
    data=aref$plot.references)+
  ggplot2::scale_x_log10(
    "N",
    breaks=meas[, 10^seq(
      ceiling(min(log10(N))),
      floor(max(log10(N))))])+
  ggplot2::scale_y_log10(
    "median line, min/max band")+
  ggplot2::geom_point(ggplot2::aes(
    N, unit.value, color=expr.name),
    data=pred,
    shape=21,
    fill="white")+
  directlabels::geom_dl(ggplot2::aes(
    N, unit.value, 
    label=label,
    color=expr.name),
    data=pred,
    method="top.polygons")+
  directlabels::geom_dl(ggplot2::aes(
    N, reference, 
    label=fun.name),
    data=aref$plot.references,
    method="bottom.polygons")+
  ggplot2::theme(legend.position="none")+
  geom_blank(aes(N, sec), data=data.frame(N=1000, sec=6, unit="seconds"))
png("figure-splits-time.png", width=8, height=4.2, units="in", res=200)
print(gg)
dev.off()

gg <- ggplot2::ggplot()+
  ggplot2::theme_bw()+
  ggplot2::facet_grid(unit ~ ., scales="free")+
  ggplot2::geom_hline(ggplot2::aes(
    yintercept=unit.value),
    color="grey",
    data=one)+
  ggplot2::geom_text(ggplot2::aes(
    100, unit.value, label=paste0(unit,"=",unit.value)),
    color="grey50",
    hjust=0,
    vjust=1.2, 
    data=one)+
  ggplot2::geom_ribbon(ggplot2::aes(
    N, ymin=min, ymax=max, fill=expr.name),
    data=meas[unit=="seconds"],
    alpha=0.5)+
  ggplot2::geom_line(ggplot2::aes(
    N, empirical, color=expr.name),
    size=2,
    data=meas)+
  geom_line(aes(
    N, reference, group=expr.name),
    data=aref$plot.references)+
  ggplot2::scale_x_log10(
    "N = number of data to segment",
    breaks=meas[, 10^seq(
      ceiling(min(log10(N))),
      floor(max(log10(N))))])+
  ggplot2::scale_y_log10(
    "median line, min/max band\n(10 timings per data size N)")+
  ggplot2::geom_point(ggplot2::aes(
    N, unit.value, color=expr.name),
    data=pred,
    shape=21,
    fill="white")+
  directlabels::geom_dl(ggplot2::aes(
    N, unit.value, 
    label=label,
    color=expr.name),
    data=pred,
    method="top.polygons")+
  directlabels::geom_dl(ggplot2::aes(
    N, reference, 
    label=sprintf("$O(%s)$", fun.latex)),
    data=aref$plot.references,
    ## method=directlabels::polygon.method("bottom", custom.colors=list(
    ##   text.color="black",
    ##   colour="white"))
    method="bottom.polygons", color="white"
    )+
  ggplot2::theme(legend.position="none")+
  geom_blank(aes(N, sec), data=data.frame(N=1000, sec=6, unit="seconds"))
tikzDevice::tikz("figure-splits-time.tex", width=6, height=3.2, standAlone=TRUE)
print(gg)
dev.off()
system("pdflatex figure-splits-time")
system("evince figure-splits-time.pdf &")
