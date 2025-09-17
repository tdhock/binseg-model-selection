library(animint2)
library(data.table)
n.segments <- 5
set.seed(1)
(seg.mean.vec <- rnorm(n.segments))
data.per.segment <- 10
data.mean.vec <- rep(seg.mean.vec, each=data.per.segment)
n.data <- length(data.mean.vec)
total.validation.sets <- 100
n.folds.vec <- 2
prop.valid.vec <- 1/n.folds.vec
set_data_dt_list <- list()
residual_dt_list <- list()
loss_dt_list <- list()
selection_dt_list <- list()
selection_error_list <- list()
for(data_seed in 1){
  set.seed(data_seed)
  data.vec <- rnorm(n.data, data.mean.vec, 0.1)
  split_data <- data.table(
    data_seed, signal=data.vec,
    position=as.numeric(seq_along(data.vec))
  )[, pos := position]
  setkey(split_data, position, pos)
  is.valid.vec.list <- list()
  for(n.folds in n.folds.vec){
    uniq.folds <- 1:n.folds
    n.seeds <- total.validation.sets/n.folds
    split_type <- sprintf("%d-fold", n.folds)
    for(cv_split_seed in 1:n.seeds){
      set.seed(cv_split_seed)
      fold.vec <- sample(rep(uniq.folds, l=n.data))
      for(valid.fold in uniq.folds){
        is.valid.vec.list[[split_type]][[paste(cv_split_seed, valid.fold)]] <-
          fold.vec==valid.fold
      }
    }
  }
  for(prop.valid in prop.valid.vec){
    split_type <- sprintf("%d%% valid", 100*prop.valid)
    prop.vec <- c(subtrain=1-prop.valid, validation=prop.valid)
    for(split_seed in 1:total.validation.sets){
      set.seed(split_seed)
      is.valid.vec.list[[split_type]][[split_seed]] <- binsegRcpp::random_set_vec(
        n.data, prop.vec) == "validation"
    }
  }
  for(split_type in names(is.valid.vec.list)){
    for(validation_set_i in seq_len(total.validation.sets)){
      is.valid <- is.valid.vec.list[[split_type]][[validation_set_i]]
      set(split_data, j="set", value=ifelse(is.valid, "validation", "subtrain"))
      set_data_dt_list[[paste(
        split_type, validation_set_i
      )]] <- data.table(
        split_type, validation_set_i,
        split_data)
      bs.model <- split_data[, binsegRcpp::binseg_normal(
        signal, is.validation.vec=is.valid)]
      coef_dt <- coef(bs.model, 1:nrow(bs.model$splits))
      setkey(coef_dt, start.pos, end.pos)
      over_dt <- foverlaps(coef_dt, split_data)
      residual_dt_list[[paste(
        data_seed, split_type, validation_set_i
      )]] <- over_dt[, .(
        data_seed, split_type, validation_set_i,
        position, segments, signal, mean, set
      )]
      loss_dt_list[[paste(
        data_seed, split_type, validation_set_i
      )]] <- data.table(
        data_seed, split_type, validation_set_i,
        bs.model$splits[, rbind(
          data.table(set="subtrain", segments, loss),
          data.table(set="validation", segments, loss=validation.loss)
        )]
      )[, point := ifelse(loss==min(loss), "min", "other"), by=set][]
      selection_dt_list[[paste(
        data_seed, split_type, validation_set_i
      )]] <- data.table(
        data_seed, split_type, validation_set_i,
        segments=bs.model$splits[which.min(validation.loss), segments]
      )
    }
  }
}

(selection_dt <- rbindlist(selection_dt_list)[
, group := rep(rep(1:n.seeds, each=n.folds), l=.N)
][])
N_order_seeds <- 100
selection_wide_list <- list()
for(order_seed in seq_len(N_order_seeds)){
  set.seed(order_seed)
  group_ord <- data.table(group=sample(n.seeds))[, split_i := 1:.N]
  selection_ord <- selection_dt[group_ord, on="group"]
  selection_cum <- data.table(num_pairs=seq_len(n.seeds))[, {
    selection_ord[split_i<=num_pairs, .(count=.N), by=.(split_type,segments)]
  }, by=num_pairs][
  , percent := 100*count/sum(count), by=.(num_pairs,split_type)][
  , tooltip := sprintf(
  "%d segments selected in %d/%d=%.1f%% validation sets, split type %s",
  segments, count, num_pairs, percent, split_type)
  ][]
  (selection_most_freq <- selection_cum[, .SD[count==max(count)], by=.(num_pairs,split_type)])
  selection_most_freq[segments==n.segments]
  selection_wide_list[[order_seed]] <- data.table(order_seed, dcast(
    selection_most_freq,
    split_type + num_pairs ~ .,
    list(min, max),
    value.var="segments"))
}

(selection_wide <- rbindlist(selection_wide_list))
normalize <- function(x)(x-min(x))/(max(x)-min(x))
log.pos.neg <- function(x)sign(x)*normalize(log10(abs(x)))
selection_wide[segments_min!=n.segments, let(
  segments_err = log.pos.neg(segments_min-n.segments)
)][, let(
  segments_range = ifelse(segments_min<segments_max, paste0(segments_min,"-",segments_max), segments_min),
  first_good_pair={
    maybe_empty <- num_pairs[segments_min!=n.segments]
    1+if(length(maybe_empty))max(maybe_empty) else 0
}), by=.(order_seed, split_type)][]
setkey(selection_wide, split_type, first_good_pair, order_seed, num_pairs)
selection_wide[, let(
  good_order = as.integer(factor(order_seed, unique(order_seed)))
), by=split_type]
good_dt <- unique(selection_wide[, .(split_type, order_seed, first_good_pair)])
good_stats <- dcast(
  good_dt,
  split_type ~ .,
  list(mean, sd),
  value.var="first_good_pair")
gg <- ggplot()+
  ggtitle(sprintf(
    "Simulation: %d data, %d true segments", n.data, n.segments))+
  geom_tile(aes(
    good_order, num_pairs, fill=segments_err),
    color=NA,
    data=selection_wide)+
  geom_label(aes(
    50, 50, label=sprintf(
      "%.1f±%.1f pairs of splits required to select true number of segments", first_good_pair_mean, first_good_pair_sd)),
    data=good_stats)+
  scale_fill_gradient2("Normalized\ndifference\nbetween\ntrue\nand\nselected\nnumber\nof\nsegments")+
  facet_grid(split_type ~ ., labeller=label_both)+
  scale_x_continuous(
    sprintf(
      "Random ordering of %d subtrain/validation splits",
      total.validation.sets),
    breaks=seq(0,100,by=10))+
  scale_y_continuous(
    "Pairs of splits used to select number of segments")+
  coord_equal()
png("figure-most-frequently-selected-orderings.png", width=7, height=6, units="in", res=200)
print(gg)
dev.off()
