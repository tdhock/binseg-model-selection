ares <- atime::atime(
  N=2^seq(2, 20),
  setup={
    x <- 1:N
  },
  verbose=TRUE,
  ## binsegRcpp={
  ##   L <- binsegRcpp::binseg_normal(x)
  ##   data.frame(splits=with(L$splits, sum(c(before.size, after.size)-1, na.rm=TRUE)))
  ## },
  #result=function(L)with(data.frame(L$res), data.frame(splits=sum(e-s))),
  sbs=wbs::sbs(x),
  seconds.limit=1,
  ## wbs "…does not lead to a significant increase in computational complexity" https://arxiv.org/abs/1411.0858
  "wbs(M=5000)"=wbs::wbs(x, 5000),
  "wbs(M=10)"=wbs::wbs(x, 10),
  "wbs(M=N/2)"=wbs::wbs(x, N/2),
  "wbs(M=N/10)"=wbs::wbs(x, N/10))
plot(ares)

saveRDS(ares, "figure-wbs-sim-data.rds")


