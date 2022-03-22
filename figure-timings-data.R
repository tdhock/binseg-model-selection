library(data.table)
timing.dt.list <- list()
for(N.data.exp in 2:14){
  N.data <- 2^N.data.exp
  max.segs <- N.data/2
  max.changes <- max.segs-1
  print(N.data)
  data.list <- list(
    best=1:N.data,
    worst=rep(c(0,1),l=N.data))
  for(case in names(data.list)){
    data.vec <- data.list[[case]]
    N.df <- microbenchmark::microbenchmark(changepoint={
      cpt.fit <- changepoint::cpt.mean(data.vec, method="BinSeg", Q=max.changes)
    }, binsegRcpp={
      binseg.fit <- binsegRcpp::binseg_normal(data.vec, max.segs)
    }, times=5)
    timing.dt.list[[paste(N.data, case)]] <- data.table(
      N.data, case, N.df)
  }
}
timing.dt <- do.call(rbind, timing.dt.list)
data.table::fwrite(timing.dt, "figure-timings-data.csv")
