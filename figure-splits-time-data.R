fun.list <- list(
  best=function(N)1:N,
  worst=function(N)rep(0:1, length.out=N))
expr.list <- atime::atime_grid(
  list(case=names(fun.list)),
  BinSeg={
    x <- data.list[[case]]
    L <- wbs::sbs(x)
    with(data.frame(L$res), data.frame(splits=sum(e-s)))
  })
ares <- atime::atime(
  N=2^seq(2, 20),
  setup={
    data.list <- sapply(fun.list, function(f)f(N), simplify=FALSE)
  },
  verbose=TRUE,
  expr.list=expr.list,
  seconds.limit=1)
plot(ares)

saveRDS(ares, "figure-splits-time-data.rds")


