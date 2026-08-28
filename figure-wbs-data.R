ares <- atime::atime(
  N=2^seq(2, 20),
  setup={
    x <- rep(0:1, length.out = N)
  },
  verbose=TRUE,
  ## binsegRcpp={
  ##   L <- binsegRcpp::binseg_normal(x)
  ##   data.frame(splits=with(L$splits, sum(c(before.size, after.size)-1, na.rm=TRUE)))
  ## },
  sbs={
    L <- wbs::sbs(x)
    with(data.frame(L$res), data.frame(splits=sum(e-s)))
  },
  wbs={
    L <- wbs::wbs(x, N/2)
    with(data.frame(L$res), data.frame(splits=sum(e-s)))
  })

data(neuroblastoma, package="neuroblastoma")
library(data.table)
k <- c("profile.id","chromosome")
nb.dt <- data.table(neuroblastoma$profiles, key=k)
count.dt <- nb.dt[, .(rows=.N), by=k]
wbs.min.N <- 4 #wbs errors for smaller N.
target.dt <- count.dt[, .(
  target=as.integer(10^seq(
    log10(wbs.min.N),
    log10(max(rows)),
    length.out = 50))
)]
#target.dt <- data.table(target=2)
(join.dt <- count.dt[target.dt, on=.(rows=target), roll="nearest", mult="first"])
(uniq.dt <- count.dt[unique(join.dt[, .(profile.id, chromosome)]), on=k])
# 48 sequences from 4 to 5937 data.

ares <- atime::atime(
  N=uniq.dt$rows,
  setup={
    meta.dt <- uniq.dt[N==rows]
    x <- nb.dt[meta.dt, logratio, on=k]
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

saveRDS(ares, "figure-wbs-data.rds")


