source("packages.R")
Sys.setenv(RETICULATE_PYTHON=if(.Platform$OS.type=="unix")
  "/home/tdhock/.local/share/r-miniconda/envs/cs570s22/bin/python"
  else "~/Miniconda3/envs/cs570s22/python.exe")
reticulate::use_condaenv("cs570s22", required=TRUE)
ruptures <- reticulate::import("ruptures")
x <- c(0,0.3,0.2,0.1, 10,11,12,13)
N.data <- length(x)
x.mat <- matrix(x)
n_bkps <- 3
changepoint::cpt.meanvar(x, method="BinSeg", Q=n_bkps)
bs.fit <- binsegRcpp::binseg("meanvar_norm", x, n_bkps+1)
seg.dt <- coef(bs.fit, 4L)
seg.dt[, n.data := end-start+1]
seg.dt[, data.table(
  data=x,
  mean=dput(rep(mean,n.data)),
  var=dput(rep(var,n.data)))]
binseg_instance <- ruptures$Binseg(model="normal", jump=1L)$fit(x.mat)
bs.fit$splits$loss
## Converting ruptures loss to normal negative log likelihood.
convert <- function(rcost)(rcost+N.data*(1+log(2*pi)))/2
(b0 <- binseg_instance$predict(n_bkps=0))
myvar <- function(y)mean((y-mean(y))^2)
v <- myvar(x)
eps <- 1e-6
log(v+eps)*N.data
binseg_instance$cost$sum_of_costs(list(b0))
(b1 <- binseg_instance$predict(n_bkps=1))
convert(binseg_instance$cost$sum_of_costs(b1))
binseg_instance$predict(n_bkps=2)
binseg_instance$predict(n_bkps=3)
