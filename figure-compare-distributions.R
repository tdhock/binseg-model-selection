source("packages.R")
dput(Sys.glob("*.csv"))
timings.csv.vec <- c(
  "normal\nmean"="figure-timings-data.csv", 
  "Laplace\nmedian"="figure-timings-laplace-data.csv", 
  "normal\nmean+\nvariance"="figure-timings-meanvar_norm-data.csv", 
  poisson="figure-timings-poisson-data.csv")
timings.dt.list <- list()
for(distribution in names(timings.csv.vec)){
  timings.csv <- timings.csv.vec[[distribution]]
  timings.dt.list[[distribution]] <- data.table(
    distribution, 
    data.table::fread(timings.csv)[grepl("binsegRcpp", package)])
}
(timings.dt <- do.call(rbind, timings.dt.list))
timings.dt[, container := sub("binsegRcpp.", "", package)]

ggplot()+
  facet_grid(container ~ case, labeller=label_both)+
  geom_ribbon(aes(
    N.data, ymin=seconds_min, ymax=seconds_max, fill=distribution),
    data=timings.dt,
    alpha=0.5)+
  geom_line(aes(
    N.data, seconds_median, color=distribution),
    data = timings.dt)+
  scale_x_log10()+
  scale_y_log10()

only.multiset <- timings.dt[container=="multiset"]
ref.dt <- rbind(
  data.table(seconds=1, unit="1 second"),
  data.table(seconds=60, unit="1 minute"))
gg <- ggplot()+
  ggtitle("Comparing distributions using binsegRcpp multiset")+
  theme_bw()+
  theme(panel.spacing=grid::unit(0, "lines"))+
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
  facet_grid(. ~ case, labeller=label_both)+
  geom_ribbon(aes(
    N.data, ymin=seconds_min, ymax=seconds_max, fill=distribution),
    data=only.multiset,
    alpha=0.5)+
  geom_line(aes(
    N.data, seconds_median, color=distribution),
    data = only.multiset)+
  scale_x_log10(
    "Number of data points to segment (log scale)",
    breaks=10^seq(1, 6, by=1))+
  coord_cartesian(
    expand=FALSE,
    ylim=c(1e-4, 1e3),
    xlim=c(8, 9e7))+
  scale_y_log10(
    "Computation time (seconds, log scale)\nMedian line and min/max band over 5 timings",
    breaks=10^seq(-10,10))
dl <- directlabels::direct.label(gg, list(cex=0.6, "right.polygons"))
png("figure-compare-distributions.png", width=9, height=3.5, units="in", res=200)
print(dl)
dev.off()
