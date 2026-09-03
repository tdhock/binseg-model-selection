library(data.table)
counts.RData.vec <- Sys.glob("data/*/*/counts.RData")
iterations.dt.list <- list()
for(file.i in seq_along(counts.RData.vec)){
  counts.RData <- counts.RData.vec[[file.i]]
  complexity.csv <- sub("counts.RData", "complexity.csv", counts.RData)
  if(file.exists(complexity.csv)){
    complexity.dt <- data.table::fread(complexity.csv)
  }else{
    cat(sprintf(
      "%4d / %4d files %s\n",
      file.i, length(counts.RData.vec), counts.RData))
    (objs <- load(counts.RData))
    count.dt <- data.table(counts)
    count.dt[, bases := chromEnd-chromStart]
    complexity.dt <- count.dt[, {
      fit <- binsegRcpp::binseg(
        "poisson", coverage, weight.vec=bases)#, max.segments=19L)
      comp.dt <- binsegRcpp::get_complexity_empirical(fit$splits)
      data.table(file.i, N.data=.N, comp.dt)
    }, by=.(cell.type, sample.id)]
    data.table::fwrite(complexity.dt, complexity.csv)
  }
  complexity.dt[, `:=`(
    sum.splits=cumsum(splits),
    max.depth=cummax(depth)
  ), by=.(cell.type, sample.id, file.i, N.data)]
  complexity.dt[segments==19, max.segments := "19"]
  complexity.dt[segments==N.data, max.segments := "N.data"]
  iterations.dt.list[[file.i]] <- complexity.dt[
    !is.na(max.segments), .(
      cell.type, sample.id, file.i, N.data,
      max.segments, sum.splits, max.depth)]
}
iterations.dt <- do.call(rbind, iterations.dt.list)
data.table::fwrite(iterations.dt, "figure-mcgill-iterations-data.csv")

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
fwrite(bound.dt, "figure-mcgill-iterations-bound.csv")
