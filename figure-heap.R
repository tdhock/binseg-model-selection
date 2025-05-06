(objs <- load("figure-heap-data.RData"))
requireNamespace("atime")
atime_list$measurements[, expr.name := sub("container=", "", expr.name)]
plot(atime_list)

fun_list <- atime:::references_funs[c("N","N \\log N", "N^2", "N^3")]
fun_list[["N^1.5"]] <- function(N)1.5 * log10(N)
ref_list <- atime::references_best(atime_list, fun_list)
ref_list$plot.references <- ref_list$references[(
  fun.name %in% c("N","N^2") & unit=="kilobytes"
)|(
  fun.name != "N" & unit=="seconds"
)]
##ref_list$measurements <- ref_list
png("figure-heap-refs-binsegRcpp.png",width=6,height=3.5,units="in",res=200)
plot_algos <- function(regex){
  regex_list <- ref_list
  regex_list$measurements <- ref_list$measurements[grepl(regex,expr.name)]
  regex_list$plot.references <- ref_list$plot.references[grepl(regex,expr.name)]
  plot(regex_list)+
    geom_blank(aes(
      N, empirical),
      data=ref_list$measurements[,.(N=1000,empirical,unit)])+
    scale_x_log10(
      breaks=10^seq(1,7,by=1))+
    theme(
      axis.text.x=element_text(angle=30,hjust=1))
}
plot_algos("list|multiset")
dev.off()

png("figure-heap-refs-other.png",width=6,height=3.5,units="in",res=200)
plot_algos("change|rupt")
dev.off()

tikz_algos <- function(regex){
  regex_list <- ref_list
  regex_list$measurements <- ref_list$measurements[grepl(regex,expr.name)]
  regex_list$plot.references <- ref_list$plot.references[grepl(regex,expr.name)][
   ,fun.name := sprintf("$%s$", sub("1.5","{1.5}",fun.latex))][]
  plot(regex_list)+
    geom_blank(aes(
      N, empirical),
      data=ref_list$measurements[,.(N=1000,empirical,unit)])+
    scale_x_log10(
      breaks=10^seq(1,7,by=1))+
    theme(
      axis.text.x=element_text(angle=30,hjust=1))
}
library(tikzDevice)
tikz("figure-heap-refs-binsegRcpp.tex",width=6,height=3)
tikz_algos("list|multiset")
dev.off()

tikz("figure-heap-refs-other.tex")
tikz_algos("change|rupt")
dev.off()

pred_list <- predict(ref_list)
png("figure-heap-pred.png",width=8,height=3,units="in",res=200)
pkg.colors <- c(
  changepoint="#E41A1C",
  ruptures="#377EB8",
  "binsegRcpp\nmultiset"="black",
  "binsegRcpp\npriority_queue"="grey30",
  "binsegRcpp\nlist"="grey50")
plot(pred_list)+
  facet_null()+
  scale_y_log10(
    "Computation time (seconds, log scale)",
    breaks=10^seq(-3,1),
    limits=c(NA,1e2))+
  scale_x_log10(
    "N = number of data to segment (max segments = N/2)",
    breaks=10^seq(1,7),
    limits=c(NA, 1e8))+
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
