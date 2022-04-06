source("packages.R")
timing.stats <- data.table::fread("figure-timings-meanvar_norm-data.csv")
disp.pkg <- c(
  changepoint="changepoint.array",
  ruptures="ruptures.LRU cache",
  "wbs::sbs"="wbs::sbs.recursion",
  "fpop::multiBinSeg"="fpop::multiBinSeg.heap")
timing.stats[, new.pkg := ifelse(
  package %in% names(disp.pkg), disp.pkg[package], package)]
timing.stats[, Package := sub("[.]", "\n", new.pkg, perl=TRUE)]
##timing.stats[, Package := sub("[.]|(?<=fpop)::", "\n", package, perl=TRUE)]
(total.minutes <- timing.stats[, sum(seconds_median*seconds_timings) / 60])
timing.stats[, min.loss := min(loss), by=.(case,N.data)]
timing.stats[loss==min.loss, .(N.data, max.segs, case, Package, loss, min.loss)]
timing.stats[loss>min.loss, .(N.data, max.segs, case, Package, loss, min.loss)]
ref.dt <- rbind(
  data.table(seconds=1, unit="1 second"),
  data.table(seconds=60, unit="1 minute"))
gg <- ggplot()+
  theme_bw()+
  theme(panel.spacing=grid::unit(0, "lines"))+
  facet_grid(. ~ case, labeller=label_both)+
  geom_line(aes(
    N.data, loss/N.data, color=Package),
    data=timing.stats)+
  coord_cartesian(
    xlim=c(NA, 1e8))+
  scale_x_log10(
    "Number of data points to segment (log scale)",
    breaks=10^seq(1, 6, by=1))+
  scale_y_log10(
    "Mean squared error (log scale)")
(dl <- directlabels::direct.label(gg, list(cex=0.6, "last.polygons")))
png("figure-timings-meanvar_norm-loss.png", width=9, height=3.5, units="in", res=200)
print(dl)
dev.off()

gg <- ggplot()+
  ggtitle("Normal change in mean and variance")+
  theme_bw()+
  theme(
    legend.position="none",
    panel.spacing=grid::unit(0, "lines"))+
  facet_grid(. ~ case, labeller=label_both)+
  geom_hline(aes(
    yintercept=seconds),
    data=ref.dt,
    color="grey")+
  geom_text(aes(
    10, seconds, label=unit),
    data=ref.dt,
    size=3,
    hjust=0,
    vjust=1.1,
    color="grey50")+
  directlabels::geom_dl(aes(
    N.data, seconds_median, color=Package, label=Package),
    data=timing.stats,
    method=list(cex=0.6, "last.polygons"))+
  geom_ribbon(aes(
    N.data, ymin=seconds_min, ymax=seconds_max, fill=Package),
    alpha=0.5,
    data=timing.stats)+
  geom_line(aes(
    N.data, seconds_median, color=Package),
    data=timing.stats)+
  scale_x_log10(
    "Number of data points to segment (log scale)",
    breaks=10^seq(1, 6, by=1))+
  coord_cartesian(
    expand=FALSE,
    ylim=c(1e-4, 1e3),
    xlim=c(8, 3e7))+
  scale_y_log10(
    "Computation time (seconds, log scale)\nMedian line and min/max band over 5 timings",
    breaks=10^seq(-10,10))
png("figure-timings-meanvar_norm.png", width=9, height=3.5, units="in", res=200)
print(gg)
dev.off()
