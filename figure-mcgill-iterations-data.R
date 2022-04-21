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
