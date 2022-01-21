library(ggplot2)
library(data.table)
future::plan("multisession")
OneExperiment <- function(experiment.i){
  experiment.row <- experiment.dt[experiment.i]
  seg.mean.vec <- 1:experiment.row$n.segments
  data.per.segment <- 10
  data.mean.vec <- rep(seg.mean.vec, each=data.per.segment)
  n.data <- length(data.mean.vec)
  n.validation.sets <- 100
  n.folds.vec <- c(10, 2)
  prop.valid.vec <- 1/n.folds.vec
  n.exp.dt <- data.table(n.exp=c(-1, -0.5, 0))
  set.seed(experiment.row$data.seed)
  data.vec <- rnorm(n.data, data.mean.vec, 0.1)
  full.model <- binsegRcpp::binseg_normal(data.vec)
  full.selection <- n.exp.dt[, {
    norm.dt <- full.model$splits[, data.table(
      norm.loss=loss*n.data^n.exp, segments)]
    penaltyLearning::modelSelection(
      norm.dt, "norm.loss", "segments")
  }, by=n.exp]
  is.valid.vec.list <- list()
  for(n.folds in n.folds.vec){
    uniq.folds <- 1:n.folds
    n.seeds <- n.validation.sets/n.folds
    split.type <- sprintf("%d-fold %d times", n.folds, n.seeds)
    for(seed in 1:n.seeds){
      set.seed(seed)
      fold.vec <- sample(rep(uniq.folds, l=n.data))
      for(valid.fold in uniq.folds){
        is.valid.vec.list[[split.type]][[paste(seed, valid.fold)]] <-
          fold.vec==valid.fold
      }
    }
  }
  for(prop.valid in prop.valid.vec){
    split.type <- sprintf("%d%% %d times", 100*prop.valid, n.validation.sets)
    prop.vec <- c(subtrain=1-prop.valid, validation=prop.valid)
    for(split.i in 1:n.validation.sets){
      set.seed(split.i)
      is.valid.vec.list[[split.type]][[split.i]] <- binsegRcpp::random_set_vec(
        n.data, prop.vec) == "validation"
    }
  }
  loss.dt <- CJ(split.i=1:n.validation.sets, type=names(is.valid.vec.list))[, {
    is.valid <- is.valid.vec.list[[type]][[split.i]]
    n.train <- n.data-sum(is.valid)
    bs.model <- binsegRcpp::binseg_normal(data.vec, is.validation.vec=is.valid)
    bs.model$splits[, data.table(
      segments,
      loss,
      n.train,
      validation.loss)]
  }, by=.(split.i, type)]
  loss.exp.dt <- n.exp.dt[, loss.dt[, .(
    split.i, type,
    segments,
    norm.loss=loss*n.train^n.exp,
    validation.loss
  )], by=n.exp]
  selection.dt <- loss.exp.dt[, penaltyLearning::modelSelection(
    .SD, "norm.loss", "segments"),
    by=.(n.exp, split.i, type)]
  selection.dt[, mid.log.lambda := ifelse(
    max.log.lambda==Inf, min.log.lambda+1, ifelse(
      min.log.lambda==-Inf, max.log.lambda-1,
      (min.log.lambda+max.log.lambda)/2))]
  selection.dt[
  , is.min := validation.loss==min(validation.loss), by=.(n.exp, split.i, type)]
  mean.selection <- selection.dt[, {
    mid.dt <- data.table(log.lambda=unique(sort(mid.log.lambda)))
    .SD[mid.dt, .(
      log.penalty=log.lambda,
      splits=.N,
      validation.loss=mean(validation.loss),
      prop.min=mean(is.min)
    ), by=.EACHI,
    on=.(min.log.lambda < log.lambda, max.log.lambda > log.lambda)]
  }, by=.(n.exp, type)]
  vline.dt <- mean.selection[
  , .SD[validation.loss == min(validation.loss)][ceiling(.N/2)],
    by=.(n.exp, type)]
  ggplot()+
    geom_vline(aes(
      xintercept=log.penalty),
      data=vline.dt)+
    geom_point(aes(
      log.penalty, validation.loss),
      shape=1,
      data=mean.selection)+
    facet_grid(n.exp ~ type)+
    scale_y_log10()
  vline.max <- mean.selection[
  , .SD[prop.min == max(prop.min)][ceiling(.N/2)],
    by=.(n.exp, type)]
  ggplot()+
    geom_vline(aes(
      xintercept=log.penalty),
      data=vline.max)+
    geom_point(aes(
      log.penalty, prop.min),
      shape=1,
      data=mean.selection)+
    facet_grid(type ~ n.exp)
  loss.stats <- loss.dt[, .(
    mean.valid.loss=mean(validation.loss)
  ), by=.(type, segments)]
  select.each.split <- loss.dt[
  , .SD[which.min(validation.loss)],
    by=.(type, split.i)]
  selected.times <- select.each.split[, .(
    times=.N
  ), by=.(type, segments)]
  lambda.dt <- rbind(
    data.table(vline.dt, what="min loss"),
    data.table(vline.max, what="max freq")
  )[, {
    full.selection[.SD, .(
      n.exp, type, segments
    ), on=.(
      n.exp, min.log.lambda < log.penalty, max.log.lambda > log.penalty
    )]
  }, by=what]
  selected.segments <- rbind(
    lambda.dt[, .(
      method=paste(type, "penaltyExp", n.exp, what),
      selected=segments
    )],
    select.each.split[, .(
      selected=min(segments)
    ), by=.(method=paste(type, "min err, min segs"))],
    selected.times[, .(
      selected=segments[which.max(times)]
    ), by=.(method=paste(type, "min err, max times"))],
    loss.stats[, .(
      selected=segments[which.min(mean.valid.loss)]
    ), by=.(method=paste(type, "mean err, min err"))]
  )
  data.table(experiment.row, selected.segments)
}
experiment.dt <- CJ(
  data.seed=1:1000,
  n.segments=c(1, 2, 5, 10, 15, 20, 50, 100))
sim.result.list <- future.apply::future_lapply(
  1:nrow(experiment.dt), OneExperiment, future.seed=NULL)
sim.result <- do.call(rbind, sim.result.list)
data.table::fwrite(sim.result, "figure-select-segments-data.csv")
