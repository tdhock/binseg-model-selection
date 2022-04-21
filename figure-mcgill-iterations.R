library(ggplot2)
iterations.dt <- data.table::fread("figure-mcgill-iterations-data.csv")
iterations.dt[order(N.data)]

ggplot()+
  geom_point(aes(
    N.data, sum.splits, color=max.segments),
    data=iterations.dt)+
  scale_x_log10()+
  scale_y_log10()

N.data.vec <- as.integer(c(200, 10^seq(2, 5)))
bound.dt.list <- list()
for(N.data in N.data.vec){
  print(N.data)
  heuristic.dt <- data.table(
    binsegRcpp::get_complexity_best_heuristic_equal_depth_full(
      N.data=N.data, min.segment.length=1L))
  heuristic.dt[, `:=`(
    sum.splits=cumsum(splits),
    max.depth=cummax(depth)
  )]
  for(n.segments in c(19L, N.data)){
    totals <- if(N.data <= 200){
      best.worst <- binsegRcpp::get_complexity_extreme(
        N.data=N.data, min.segment.length=1L, n.segments=n.segments)
      best.worst[, .(
        sum.splits=sum(splits),
        max.depth=max(depth)
      ), by=case]
    }else data.table(
      case="worst",
      sum.splits=n.segments*(1+N.data-(n.segments+3)/2),
      max.depth=N.data-1)
    bound.dt.list[[paste(N.data, n.segments)]] <- data.table(
      N.data, n.segments, rbind(
        heuristic.dt[n.segments, data.table(
          case="best.heuristic",
          sum.splits,
          max.depth)],
        totals))
  }
}

bound.dt <- do.call(rbind, bound.dt.list)
bound.dt[, max.segments := ifelse(n.segments==19, "19", "N.data")]
size <- 1
gg <- ggplot()+
  geom_point(aes(
    N.data, sum.splits),
    shape=1,
    data=iterations.dt)+
  geom_line(aes(
    N.data, sum.splits, color=case, size=case),
    data=bound.dt)+
  facet_grid(. ~ max.segments, labeller=label_both)+
  scale_x_log10(
    "Number of data to segment")+
  scale_y_log10(
    "Number of candidate\nsplits to consider")+
  theme_bw()+
  theme(panel.spacing=grid::unit(0,"lines"))+
  coord_cartesian(ylim=c(1e2, 1e7))+
  scale_size_manual(values=c(
    best=1.5, best.heuristic=0.75, worst=0.75))
png("figure-mcgill-iterations.png", width=7, height=2, units="in", res=200)
print(gg)
dev.off()

iterations.dt[, round.N.data := 10^(round(log10(N.data)))]
iterations.dt[, num.sum.splits := as.numeric(sum.splits)]
iterations.stats <- dcast(
  iterations.dt,
  round.N.data + max.segments ~ .,
  fun.aggregate=list(median, min, max, length),
  value.var="num.sum.splits")
iterations.stats[, case := "empirical"]
case.colors <- c(
  empirical="black",
  best="#E41A1C",
  worst="#377EB8",
  best.heuristic="#4DAF4A")#, "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF")
ggplot()+
  scale_color_manual(values=case.colors)+
  scale_fill_manual(values=case.colors, guide="none")+
  geom_ribbon(aes(
    round.N.data,
    ymin=num.sum.splits_min,
    ymax=num.sum.splits_max,
    fill=case),
    alpha=0.5,
    data=iterations.stats)+
  geom_line(aes(
    round.N.data, num.sum.splits_median, color=case),
    size=size,
    data=iterations.stats)+
  facet_grid(. ~ max.segments, labeller=label_both)+
  scale_x_log10(
    "Number of data to segment")+
  scale_y_log10(
    "Number of candidate splits to consider")+
  geom_line(aes(
    N.data, sum.splits, color=case),
    size=size,
    data=bound.dt)+
  theme_bw()+
  theme(panel.spacing=grid::unit(0,"lines"))+
  coord_cartesian(ylim=c(1e2, 1e7))
