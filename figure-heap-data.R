Sys.setenv(RETICULATE_PYTHON="~/miniconda3/envs/ruptures/bin/python")
reticulate::use_condaenv("ruptures", required=TRUE)

ruptures <- reticulate::import("ruptures")
binseg_l2 <- ruptures$Binseg(min_size=1L, jump=1L)

gen_data <- function(N.data)list(
  flat=rnorm(N.data),
  best=1:N.data,
  worst=rep(c(0,1),l=N.data))
expr.list <- atime::atime_grid(list(
  container=c("priority_queue","list","multiset")),
  ##case=names(gen_data(2))),
  binsegRcpp={
    binseg.fit <- binsegRcpp::binseg(
      "mean_norm", data_vec, max.segs, container.str = container)
    sort(binseg.fit$splits$end)
  },
  expr.param.sep="\n")
atime_list <- atime::atime(
  N=2^seq(2, 30),
  ##N=2^seq(2, 4),
  setup={
    max.segs <- as.integer(N/2)
    max.changes <- max.segs-1L
    set.seed(1)
    ## data_list <- gen_data(N)
    ## mat_list <- lapply(data_list, cbind)
    data_vec <- 1:N
    data_mat <- cbind(data_vec)
  },
  ruptures={
    binseg_l2$fit(data_mat)$predict(n_bkps=max.changes)
  },
  changepoint={
    cpt.fit <- changepoint::cpt.mean(data_vec, method="BinSeg", Q=max.changes)
    sort(c(N, cpt.fit@cpts.full[max.changes,]))
  },
  verbose=TRUE,
  result=TRUE,
  seconds.limit=5,
  expr.list=expr.list)
atime_list$measurements[N<20, .(N, expr.name, paste(result))]
library(data.table)
atime_wide <- dcast(
  atime_list$measurements[!is.na(container)],
  N ~ container,
  value.var="result")
stopifnot(atime_wide[, identical(multiset, priority_queue)])
atime_list$measurements[, result := NULL][]
save(atime_list, file="figure-heap-data.RData")
