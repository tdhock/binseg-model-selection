library(ggplot2)
library(data.table)
timings.dt <- data.table::fread("figure-timings-data.csv")
timings.dt[, seconds := time/1e9]
timings.dt[, package := expr]
ggplot()+
  facet_grid(. ~ case, labeller=label_both)+
  geom_point(aes(
    N.data, seconds, color=package),
    data=timings.dt)+
  scale_x_log10()+
  scale_y_log10()

timing.stats <- data.table::dcast(
  timings.dt,
  case + package + N.data ~ .,
  list(median, min, max, trials=length),
  value.var="seconds")
ref.dt <- rbind(
  data.table(seconds=1, unit="1 second"),
  data.table(seconds=60, unit="1 minute"))
gg <- ggplot()+
  theme_bw()+
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
    N.data, ymin=seconds_min, ymax=seconds_max, fill=package),
    alpha=0.5,
    data=timing.stats)+
  geom_line(aes(
    N.data, seconds_median, color=package),
    data=timing.stats)+
  scale_x_log10(
    "Number of data points to segment",
    breaks=10^seq(1, 4, by=1))+
  coord_cartesian(
    expand=FALSE,
    ylim=c(1e-4, 1e2),
    xlim=c(8, 500000))+
  scale_y_log10(
    "Computation time (seconds)\nMedian line and min/max band over 5 timings",
    breaks=10^seq(-10,10))
dl <- directlabels::direct.label(gg, "last.polygons")
png("figure-timings.png", width=7, height=3.5, units="in", res=200)
print(dl)
dev.off()
