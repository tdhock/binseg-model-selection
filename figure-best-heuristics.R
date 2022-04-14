source("packages.R")

min.seg.len <- 2L
n.segs <- 4L
compare.dt.list <- list()
for(N.data in as.integer(10^seq(1, 3, by=0.25))){
  print(N.data)
  breadth.vec <- binsegRcpp:::get_best_heuristic_equal(
    N.data, min.seg.len)
  depth.df <- binsegRcpp:::depth_first_interface(
    N.data, min.seg.len)
  node.dt <- binsegRcpp:::get_best_optimal(
    N.data, min.seg.len, n.segs)
  first.splits <- binsegRcpp:::size_to_splits(N.data, min.seg.len)
  worst.splits <- c(
    if(1 <= first.splits)seq(first.splits, 1, by=-min.seg.len),
    0)[1:n.segs]
  best.dt <- node.dt[, .(
    candidates=sum(binsegRcpp:::size_to_splits(s, min.seg.len))
  ), by=.(parent,level)]
  N.dt <- rbind(
    data.table(algo="worst", candidates=sum(worst.splits)),
    data.table(algo="best", candidates=sum(best.dt$candidates)),
    data.table(algo="breadth", candidates=sum(breadth.vec[1:n.segs])),
    data.table(algo="depth", candidates=sum(depth.df$splits[1:n.segs])))
  compare.dt.list[[paste(N.data)]] <- data.table(
    N.data, N.dt)
}
compare.dt <- do.call(rbind, compare.dt.list)

gg <- ggplot()+
  ggtitle(sprintf(
    "Min segment length = %d, N segments = %d",
    min.seg.len, n.segs))+
  geom_line(aes(
    N.data, candidates, color=algo),
    data=compare.dt)+
  scale_x_log10()+
  scale_y_log10()+
  coord_cartesian(
    xlim=c(10,3000))
dl <- directlabels::direct.label(gg, "right.polygons")
png(
  "figure-best-heuristics-segs-constant.png", 
  width=4, height=3, units="in", res=200)
print(dl)
dev.off()

min.seg.len <- 2L
compare.dt.list <- list()
segs.factor <- 0.3
for(N.data in as.integer(10^seq(1, 2, by=0.25))){
  print(N.data)
  n.segs <- as.integer(N.data*segs.factor)
  breadth.vec <- binsegRcpp:::get_best_heuristic_equal(
    N.data, min.seg.len)
  depth.df <- binsegRcpp:::depth_first_interface(
    N.data, min.seg.len)
  node.dt <- binsegRcpp:::get_best_optimal(
    N.data, min.seg.len, n.segs)
  first.splits <- binsegRcpp:::size_to_splits(N.data, min.seg.len)
  worst.splits <- c(
    if(1 <= first.splits)seq(first.splits, 1, by=-min.seg.len),
    0)[1:n.segs]
  best.dt <- node.dt[, .(
    candidates=sum(binsegRcpp:::size_to_splits(s, min.seg.len))
  ), by=.(parent,level)]
  N.dt <- rbind(
    data.table(algo="worst", candidates=sum(worst.splits)),
    data.table(algo="best", candidates=sum(best.dt$candidates)),
    data.table(algo="breadth", candidates=sum(breadth.vec[1:n.segs])),
    data.table(algo="depth", candidates=sum(depth.df$splits[1:n.segs])))
  compare.dt.list[[paste(N.data)]] <- data.table(
    N.data, N.dt)
}
compare.dt <- do.call(rbind, compare.dt.list)
gg <- ggplot()+
  ggtitle(sprintf(
    "Min segment length = %d, N segments = %.1f*N data",
    min.seg.len, segs.factor))+
  geom_line(aes(
    N.data, candidates, color=algo),
    data=compare.dt)+
  scale_x_log10()+
  scale_y_log10()+
  coord_cartesian(
    xlim=c(10,150))
dl <- directlabels::direct.label(gg, "right.polygons")
png("figure-best-heuristics.png", width=5, height=3, units="in", res=200)
print(dl)
dev.off()
