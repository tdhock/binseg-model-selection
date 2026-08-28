library(ggplot2)
ares <- readRDS("figure-wbs-data.rds")
ares$unit.col.vec <- c(seconds="median")
png("figure-wbs-sim-nb.png", width=7, height=4, units="in", res=200)
plot(ares)+
  ggtitle("Neuroblastoma data")+
  facet_null()+
  scale_y_log10(
    "Computation time (seconds)\nmedian line, min/max band")+
  scale_x_log10(
    "N = number of data to segment",
    breaks=c(4, 10, 100, 1000, 5937),
    limits=c(4, 3e4))
dev.off()

ares <- readRDS("figure-wbs-sim-data.rds")
ares$unit.col.vec <- c(seconds="median")
aref <- atime::references_best(ares)
png("figure-wbs-sim-ref.png", width=12, height=3, units="in", res=200)
plot(aref)
dev.off()

apred <- predict(aref)
png("figure-wbs-sim-pred.png", width=7, height=4, units="in", res=200)
plot(apred)+
  ggtitle("Best case synthetic data")+
  facet_null()+
  scale_y_log10(
    "Computation time (seconds)\nmedian line, min/max band")+
  scale_x_log10(
    "N = number of data to segment")
dev.off()

