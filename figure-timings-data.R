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
os <- reticulate::import("os")
os$listdir(".")
builtins <- reticulate::import_builtins()
builtins$dir(ruptures)
ruptures[["__path__"]]
n_samples <- 1000L
dim <- 1L
sigma <- 4
n_bkps <- 4L
sim.list = ruptures$pw_constant(n_samples, dim, n_bkps, noise_std=sigma)
sim.mat <- sim.list[[1]]
sim.vec <- as.numeric(sim.mat)
algo = ruptures$Pelt(model="rbf")$fit(sim.mat)
changepoint::cpt.mean(sim.vec, method="BinSeg", Q=n_bkps)
binsegRcpp::binseg_normal(sim.vec, n_bkps+1)
algo$predict(pen=10)
binseg_instance <- ruptures$Binseg(min_size=1L, jump=1L)
binseg_instance$fit(sim.list[[1]])
binseg_instance$predict(n_bkps=4L)
timing.dt.list <- list()
done.list <- list()
for(N.data.exp in 2:20){#2^20 = 1,048,576
  N.data <- 2^N.data.exp
  max.segs <- as.integer(N.data/2)
  max.changes <- max.segs-1L
  print(N.data)
  set.seed(1)
  data.list <- list(
    flat=rnorm(N.data),
    best=1:N.data,
    worst=rep(c(0,1),l=N.data))
  for(case in names(data.list)){
    data.vec <- data.list[[case]]
    cum.data.vec <- cumsum(c(0,data.vec))
    data.mat <- matrix(data.vec)
    result.list <- list()
    m.args <- do.sub(ruptures={
      binseg_instance$fit(data.mat)$predict(n_bkps=max.changes)
    }, changepoint={
      cpt.fit <- changepoint::cpt.mean(data.vec, method="BinSeg", Q=max.changes)
      sort(c(N.data,cpt.fit@cpts.full[max.changes,]))
    }, binsegRcpp.multiset={
      binseg.fit <- binsegRcpp::binseg(
        "mean_norm",data.vec, max.segs, container.str="multiset")
      sort(binseg.fit$splits$end)
    }, wbs={
      wbs.fit <- wbs::sbs(data.vec)
      split.dt <- data.table(wbs.fit$res)[order(-min.th, scale)]
      sort(split.dt[, c(N.data, cpt)][1:max.segs])
    }, binsegRcpp.list={
      binseg.fit <- binsegRcpp::binseg(
        "mean_norm",data.vec, max.segs, container.str="list")
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
        seg.mean <- (cum.data.vec[end+1]-cum.data.vec[end-seg.size+1])/seg.size
        data.mean <- rep(seg.mean, seg.size)
        sum((data.mean-data.vec)^2)
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
data.table::fwrite(timing.dt, "figure-timings-data.csv")
