library(data.table)
data(neuroblastoma, package="neuroblastoma")
nb.dt <- data.table(neuroblastoma$profiles)
rle.dt <- nb.dt[, {
  L <- rle(logratio)
  do.call(data.table, L)
}, by=.(profile.id, chromosome)]
table(rle.dt$lengths)
nb.list <- split(rle.dt, rle.dt[, paste(profile.id, chromosome)])
num.data <- sapply(nb.list, nrow)
max.segs <- 10L
nb.big <- nb.list[num.data >= max.segs]
length(nb.big)
length(nb.list)

iterations.dt.list <- list()
for(seq.i in seq_along(nb.big)){
  count.dt <- nb.big[[seq.i]]
  meta.dt <- count.dt[1, .(profile.id, chromosome)]
  complexity.csv <- meta.dt[, paste0(file.path(
    "neuroblastoma", profile.id, chromosome), ".csv")]
  pro.dir <- dirname(complexity.csv)
  dir.create(pro.dir, showWarnings = FALSE, recursive = TRUE)
  if(file.exists(complexity.csv)){
    complexity.dt <- data.table::fread(complexity.csv)
  }else{
    cat(sprintf(
      "%4d / %4d files %s\n",
      seq.i, length(nb.big), complexity.csv))
    fit <- count.dt[, binsegRcpp::binseg("mean_norm", values, weight.vec=lengths)]
    comp.dt <- binsegRcpp::get_complexity_empirical(fit$splits)
    complexity.dt <- comp.dt[, .(splits, depth)]
    data.table::fwrite(complexity.dt, complexity.csv)
  }
  complexity.dt[, `:=`(
    sum.splits=cumsum(splits),
    max.depth=cummax(depth)
  )]
  complexity.dt[max.segs, max.segments := paste(max.segs)]
  complexity.dt[.N, max.segments := "N.data"]
  iterations.dt.list[[seq.i]] <- data.table(
    meta.dt, N.data=nrow(count.dt),
    complexity.dt[!is.na(max.segments), .(max.segments, sum.splits, max.depth)])
}
iterations.dt <- do.call(rbind, iterations.dt.list)
data.table::fwrite(iterations.dt, "figure-neuroblastoma-iterations-data.csv")

bound.dt.list <- list()

N.data.vec <- unique(sort(iterations.dt$N.data))
for(N.data in N.data.vec){
  if(!paste(N.data, N.data) %in% names(bound.dt.list)){
    print(N.data)
    heuristic.dt <- data.table(
      binsegRcpp::get_complexity_best_heuristic_equal_depth_full(
        N.data=N.data, min.segment.length=1L))
    heuristic.dt[, `:=`(
      sum.splits=cumsum(splits),
      max.depth=cummax(depth)
    )]
    for(n.segments in as.integer(c(max.segs, N.data))){
      totals <- if(N.data <= 100){
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
}
bound.dt <- do.call(rbind, bound.dt.list)
bound.dt[, max.segments := ifelse(n.segments==max.segs, max.segs, "N.data")]
data.table::fwrite(bound.dt, "figure-neuroblastoma-iterations-bounds.csv")
