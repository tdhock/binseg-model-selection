source("packages.R")
Sys.setenv(RETICULATE_PYTHON=if(.Platform$OS.type=="unix")
  "/home/tdhock/.local/share/r-miniconda/envs/cs570s22/bin/python"
  else "~/Miniconda3/envs/cs570s22/python.exe")
reticulate::use_condaenv("cs570s22", required=TRUE)
ruptures <- reticulate::import("ruptures")
x <- c(0,0.3,0.2,0.1, 10,11,12,13)
N.data <- length(x)
x.mat <- matrix(x)
max_bkps <- 3
changepoint::cpt.mean(x, method="BinSeg", Q=max_bkps)
mbs.fit <- fpop::multiBinSeg(x,max_bkps)
bs.fit <- binsegRcpp::binseg("mean_norm", x, max_bkps+1)

binseg_instance <- ruptures$Binseg(model="l2", jump=1L)$fit(x.mat)
(b0 <- binseg_instance$predict(n_bkps=0))
binseg_instance$cost$sum_of_costs(as.list(b0))
binseg_instance$predict(n_bkps=1)
binseg_instance$predict(n_bkps=2)
ruptures.loss <- sapply(seq(0,n_bkps), function(n_bkps){
  end.vec <- binseg_instance$predict(n_bkps=n_bkps)
  tryCatch({
    binseg_instance$cost$sum_of_costs(as.list(end.vec))
  }, error=function(e){
    NA_real_
  })
})
sbs.fit <- wbs::sbs(x)
cpt.fit <- changepoint::cpt.mean(x, method="BinSeg", Q=max_bkps)
changepoint::logLik(cpt.fit)
changepoint::cpts(cpt.fit)
cpt.fit@cpts.full
sapply(c("binsegRcpp", "changepoint", "wbs", "fpop"), function(p)paste(packageVersion(p)))
ruptures$version$version
rbind(
  ruptures=ruptures.loss,
  fpop=cumsum(c(sum((x-mean(x))^2), mbs.fit$J.est)),
  binsegRcpp=bs.fit$splits$loss)
