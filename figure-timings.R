source("packages.R")
timing.stats <- data.table::fread("figure-timings-data.csv")
timing.stats[, Package := sub("[.]", "\n", package)]
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
    xlim=c(10, 5e5))+
  scale_x_log10(
    "Number of data points to segment",
    breaks=10^seq(1, 4, by=1))+
  scale_y_log10(
    "Mean squared error")
(dl <- directlabels::direct.label(gg, list(cex=0.7, "last.polygons")))
png("figure-timings-loss.png", width=7, height=3.5, units="in", res=200)
print(dl)
dev.off()

gg <- ggplot()+
  theme_bw()+
  theme(panel.spacing=grid::unit(0, "lines"))+
  facet_grid(. ~ case, labeller=label_both)+
  geom_hline(aes(
    yintercept=seconds),
    data=ref.dt,
    color="grey")+
  geom_text(aes(
    10, seconds, label=unit),
    data=ref.dt,
    hjust=0,
    vjust=1.1,
    color="grey50")+
  geom_ribbon(aes(
    N.data, ymin=seconds_min, ymax=seconds_max, fill=Package),
    alpha=0.5,
    data=timing.stats)+
  geom_line(aes(
    N.data, seconds_median, color=Package),
    data=timing.stats)+
  scale_x_log10(
    "Number of data points to segment",
    breaks=10^seq(1, 4, by=1))+
  coord_cartesian(
    expand=FALSE,
    ylim=c(1e-4, 7e2),
    xlim=c(8, 8e7))+
  scale_y_log10(
    "Computation time (seconds)\nMedian line and min/max band over 5 timings",
    breaks=10^seq(-10,10))
(dl <- directlabels::direct.label(gg, list(cex=0.7, "last.polygons")))
png("figure-timings.png", width=7, height=3.5, units="in", res=200)
print(dl)
dev.off()
