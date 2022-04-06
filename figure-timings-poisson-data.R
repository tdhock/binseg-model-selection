source("packages.R")

seconds.limit <- 100
do.sub <- function(...){
  mcall <- match.call()
  L <- as.list(mcall[-1])
  for(arg.name in names(L)){
    maybe.lang <- L[[arg.name]]
    if(is.language(maybe.lang)){
      L[[arg.name]] <- substitute(
        result.list[[NAME]] <- EXPR,
        list(NAME=arg.name, EXPR=maybe.lang))
    }
  }
  L
}

timing.dt.list <- list()
done.list <- list()
for(N.data.exp in 2:20){#2^20 = 1,048,576
  N.data <- 2^N.data.exp
  max.segs <- as.integer(N.data/2)
  max.changes <- max.segs-1L
  print(N.data)
  set.seed(1)
  data.list <- list(
    flat=rpois(N.data,100),
    best=1:N.data,
    worst=rep(c(0,1),l=N.data))
  for(case in names(data.list)){
    data.vec <- data.list[[case]]
    cum.data.vec <- cumsum(c(0,data.vec))
    data.mat <- matrix(data.vec)
    result.list <- list()
    m.args <- do.sub("changepoint"={
      cpt.fit <- changepoint::cpt.meanvar(
        data.vec, method="BinSeg", Q=max.changes, test.stat="Poisson")
      sort(c(N.data,cpt.fit@cpts.full[max.changes,]))
    }, blockcpd={
      bcpd.fit <- blockcpd::fit_blockcpd(
        t(data.mat), method="hierseg", family="poisson", lambda=0)
      c(bcpd.fit$changepoints, N.data)
    }, "binsegRcpp.multiset"={
      binseg.fit <- binsegRcpp::binseg(
        "poisson",data.vec, max.segs, container.str="multiset")
      sort(binseg.fit$splits$end)
    }, "binsegRcpp.list"={
      binseg.fit <- binsegRcpp::binseg(
        "poisson",data.vec, max.segs, container.str="list")
      sort(binseg.fit$splits$end)
    }, times=5)
    m.args[names(done.list[[case]])] <- NULL
    if(length(m.args) > 1){
      N.df <- do.call(microbenchmark::microbenchmark, m.args)
      N.dt <- data.table(N.df)
      N.dt[, seconds := time/1e9]
      N.dt[, package := paste(expr)]
      N.stats <- data.table::dcast(
        N.dt,
        package ~ .,
        list(median, min, max, timings=length),
        value.var="seconds")
      N.stats[, computed.ends := result.list[package] ]
      N.stats[, valid.ends := lapply(computed.ends, function(ends){
        u <- unique(ends)
        u[u > 0 & u <= N.data]
      })]
      loss <- sapply(result.list, function(end){
        seg.size <- diff(c(0,end))
        seg.sum <- cum.data.vec[end+1]-cum.data.vec[end-seg.size+1]
        sum(seg.sum-log(seg.sum/seg.size))
      })
      N.stats[, loss := loss[package] ]
      done.pkgs <- N.stats[seconds_median > seconds.limit, paste(package)]
      done.list[[case]][done.pkgs] <- TRUE
      timing.dt.list[[paste(N.data, case)]] <- data.table(
        N.data, max.segs, case, N.stats)
    }
  }
}

timing.dt <- do.call(rbind, timing.dt.list)
timing.dt[N.data > 100, `:=`(computed.ends=NA, valid.ends=NA)]
data.table::fwrite(timing.dt, "figure-timings-poisson-data.csv")
