(objs <- load("figure-heap-data.RData"))
requireNamespace("atime")
atime_list$measurements[, expr.name := sub("container=", "", expr.name)]
plot(atime_list)

fun_list <- atime:::references_funs[c("N","N \\log N", "N^2", "N^3")]
myexp <- 1.4
fun_list[[paste0("N^",myexp)]] <- function(N)myexp * log10(N)
ref_list <- atime::references_best(atime_list, fun_list)
ref_list$plot.references <- ref_list$references[(
  fun.name %in% c("N","N^2") & unit=="kilobytes"
)|(unit=="seconds" & ( (
  grepl("binsegRcpp", expr.name) & fun.name %in% c("N log N", "N^2")
)|(
  expr.name=="ruptures" & !fun.name%in%c("N^3","N")
)|(
  expr.name=="changepoint" & fun.name%in%c("N^2","N^3")
)))]
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
    limits=c(NA,200))+
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

tikz_algos <- function(f.tex,regex){
  trans <- function(DT)DT[,expr.name:=sub("_","\\_",expr.name,fixed=TRUE)]
  meas <- trans(ref_list$measurements[grepl(regex,expr.name)])
  ref.dt <- trans(ref_list$plot.references[grepl(regex,expr.name)])[,fun.name := sprintf("$%s$", sub(paste(myexp),paste0("{",myexp,"}"),fun.latex))][]
  ##plot(regex_list)+
  emp.color <- "black"
  ref.color <- "violet"
  gg <- ggplot()+
    facet_grid(unit ~ expr.name, scales="free")+
    ggplot2::geom_ribbon(ggplot2::aes(
      N, ymin = min, 
      ymax = max, group = expr.name),
      data = meas[unit == "seconds"],
      fill = emp.color,
      alpha = 0.5) +
    ggplot2::geom_line(ggplot2::aes(
      N, empirical, group = expr.name),
      size = 2,
      color = emp.color, 
      data = meas) +
    ggplot2::geom_line(
      ggplot2::aes(
        N, 
        reference, group = paste(expr.name, fun.name)),
      color = ref.color, 
      size = 1, data = ref.dt) +
    geom_blank(aes(
      N, empirical/5),
      data=trans(ref_list$measurements)[,.(N=1000,empirical,unit)])+
    scale_x_log10(
      "$N$ = number of data to segment",
      breaks=10^seq(1,7,by=1))+
    scale_y_log10("")+
    theme(
      axis.text.x=element_text(angle=30,hjust=1))+
    directlabels::geom_dl(ggplot2::aes(
      N, reference, 
      label = fun.name),
      data = ref.dt,
      color = ref.color, 
      method = list("bottom.polygons", directlabels::dl.trans(y=y-0.05)))
  tikz(f.tex,width=6,height=3)
  print(gg)
  dev.off()
}
library(tikzDevice)
tikz_algos("figure-heap-refs-binsegRcpp.tex", "list|multiset")
tikz_algos("figure-heap-refs-other.tex", "change|rupt")

