source("packages.R")
timing.stats <- data.table::fread("figure-timings-versions-data.csv")
timing.stats[, new.pkg := ifelse(
  package %in% names(disp.pkg), disp.pkg[package], package)]
timing.stats[, Package := sub("[.]", "\n", new.pkg, perl=TRUE)]
(total.minutes <- timing.stats[, sum(seconds_median*seconds_timings) / 60])
ref.dt <- rbind(
  data.table(seconds=1, unit="1 second"),
  data.table(seconds=60, unit="1 minute"))
gg <- ggplot()+
  ggtitle("Normal change in mean with constant variance (square loss)")+
  theme_bw()+
  theme(
    legend.position="none",
    panel.spacing=grid::unit(0, "lines"))+
  facet_grid(. ~ case, labeller=label_both)+
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
    xlim=c(8, 5e5))+
  scale_y_log10(
    "Computation time (seconds, log scale)\nMedian line and min/max band over 5 timings",
    breaks=10^seq(-10,10))
png("figure-timings-versions.png", width=9, height=3.5, units="in", res=200)
print(gg)
dev.off()
