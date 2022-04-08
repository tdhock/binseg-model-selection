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
Sys.setenv(RETICULATE_PYTHON=if(.Platform$OS.type=="unix")
  "/home/tdhock/.local/share/r-miniconda/envs/cs570s22/bin/python"
  else "~/Miniconda3/envs/cs570s22/python.exe")
reticulate::use_condaenv("cs570s22", required=TRUE)
ruptures <- reticulate::import("ruptures")
binseg_instance <- ruptures$Binseg(min_size=1L, jump=1L, model="normal")

timing.dt.list <- list()
done.list <- list()
N.data.exp <- 17
case <- "best"
for(N.data.exp in 2:20){#2^20 = 1,048,576
  N.data <- 2^N.data.exp
  max.segs <- as.integer(N.data/2)
  max.changes <- max.segs-1L
  print(N.data)
  set.seed(1)
  data.list <- list(
    flat=rnorm(N.data),
    best=1:N.data,
    worst=rep(c(0,1,10,11),l=N.data))
  for(case in names(data.list)){
    data.vec <- data.list[[case]]
    cum.data.vec <- cumsum(c(0,data.vec))
    cum.square.vec <- cumsum(c(0,data.vec^2))
    data.mat <- matrix(data.vec)
    result.list <- list()
    m.args <- do.sub("ruptures"={
      binseg_instance$fit(data.mat)$predict(n_bkps=max.changes)
    }, "changepoint"={
      cpt.fit <- changepoint::cpt.meanvar(
        data.vec, method="BinSeg", Q=max.changes)
      sort(c(N.data,cpt.fit@cpts.full[max.changes,]))
    }, blockcpd={
      bcpd.fit <- blockcpd::fit_blockcpd(
        t(data.mat), method="hierseg", family="normal", lambda=0,
        min_block_size = 2L)
      c(bcpd.fit$changepoints, N.data)
    }, "binsegRcpp.multiset"={
      binseg.fit.multiset <- binsegRcpp::binseg(
        "meanvar_norm",data.vec, max.segs, container.str="multiset")
      sort(binseg.fit.multiset$splits$end)
    }, "binsegRcpp.list"={
      binseg.fit.list <- binsegRcpp::binseg(
        "meanvar_norm",data.vec, max.segs, container.str="list")
      sort(binseg.fit.list$splits$end)
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
      nll <- function(end){
        seg.size <- diff(c(0,end))
        end.i <- end+1
        start.i <- end-seg.size+1
        seg.sum <- cum.data.vec[end.i]-cum.data.vec[start.i]
        seg.square <- cum.square.vec[end.i]-cum.square.vec[start.i]
        var.est <- (seg.square-seg.sum^2/seg.size)/seg.size
        sum((log(var.est*2*pi)+1)*seg.size/2)
      }
      loss <- sapply(result.list, nll)
      N.stats[, loss := loss[package] ]
      done.pkgs <- N.stats[seconds_median > seconds.limit, paste(package)]
      done.list[[case]][done.pkgs] <- TRUE
      timing.dt.list[[paste(N.data, case)]] <- data.table(
        N.data, max.segs, case, N.stats)
    }
  }
}

if(FALSE){
  binseg1819[["17"]] <- list(
      list=binseg.fit.list,
      multiset=binseg.fit.multiset)
  save(binseg1819,file="binseg1819.RData")
}

timing.dt <- do.call(rbind, timing.dt.list)
timing.dt[, computed.N := sapply(computed.ends, length)]
timing.dt[N.data > 100, `:=`(computed.ends=NA, valid.ends=NA)]
data.table::fwrite(timing.dt, "figure-timings-meanvar_norm-data.csv")
