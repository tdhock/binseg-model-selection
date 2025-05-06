(objs <- load("figure-heap-data.RData"))
requireNamespace("atime")
plot(atime_list)

fun_list <- atime:::references_funs[c("N","N \\log N", "N^2", "N^3")]
fun_list[["N^1.5"]] <- function(N)1.5 * log10(N)
ref_list <- atime::references_best(atime_list, fun_list)
ref_list$plot.references <- ref_list$references[(
  fun.name %in% c("N","N^2") & unit=="kilobytes"
)|(
  fun.name != "N" & unit=="seconds"
)]
png("figure-heap-refs.png",width=14,height=6,units="in",res=200)
plot(ref_list)
dev.off()

pred_list <- predict(ref_list)
png("figure-heap-pred.png",width=12,height=4.5,units="in",res=200)
pkg.colors <- c(
  changepoint="#E41A1C",
  ruptures="#377EB8",
  "binsegRcpp\ncontainer=multiset"="black",
  "binsegRcpp\ncontainer=priority_queue"="grey30",
  "binsegRcpp\ncontainer=list"="grey50")
plot(pred_list)+
  facet_null()+
  scale_y_log10(
    "Computation time (seconds, log scale)")+
  scale_x_log10(
    "N = number of data to segment (max segments = N/2)",
    breaks=10^seq(1,7),
    limits=c(NA, 1e9))+
  scale_color_manual(values=pkg.colors)+
  scale_fill_manual(values=pkg.colors)
dev.off()

norm_dt <- pred_list$prediction[, ratio := max(N)/N][, .(expr.name, ratio)][pred_list$meas, on="expr.name"][unit=="seconds"][, norm_N := N*ratio][]
library(ggplot2)
ggplot()+
  geom_line(aes(
    norm_N, empirical, color=expr.name),
    data=norm_dt)+
  scale_x_log10()+
  scale_y_log10()
