library(data.table)
n.segments <- 5
set.seed(1)
(seg.mean.vec <- rnorm(n.segments))
data.per.segment <- 20
data.mean.vec <- rep(seg.mean.vec, each=data.per.segment)
n.data <- length(data.mean.vec)
total.validation.sets <- 2
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

my_rbind <- function(L)rbindlist(L)[, Split := as.integer(factor(paste(split_type, validation_set_i)))]
loss_dt <- my_rbind(loss_dt_list)[, Total_Square_Loss := ifelse(loss<1e-10, 0, loss)][]
loss_lab <- loss_dt[segments==max(segments)]
(set_data_dt <- my_rbind(set_data_dt_list))
library(animint2)
dput(RColorBrewer::brewer.pal(Inf, "Set2"))
c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", 
  "#E5C494", "#B3B3B3")
set_colors <- c(
  subtrain="#66C2A5",
  validation="#FC8D62")
(residual_dt <- my_rbind(residual_dt_list))
## double check residuals agree with loss.
data.table(
  loss_dt[, .(segments,set,loss)],
  residual_dt[, .(
    RSS=sum((mean-signal)^2)
  ), keyby=.(set, segments)][, .(RSS)]
)[, let(diff=loss-RSS)][]
seg_path_dt <- data.table(offset=0.5*c(1,-1))[, residual_dt[, .(
  position=position-offset, mean, segments, Split, data_i=position
)], by=offset][order(segments, data_i, position)]
viz <- animint(
  data=ggplot()+
    theme_bw()+
    theme_animint(width=800)+
    scale_fill_manual(values=set_colors)+
    scale_color_manual(values=set_colors)+
    geom_segment(aes(
      position, signal, key=position,
      color=set,
      xend=position, yend=mean),
      showSelected=c("Split","segments"),
      size=1,
      help="Vertical segments show residuals (difference between data and mean)",
      data=residual_dt)+
    geom_point(aes(
      position, signal,
      fill=set,
      key=position),
      size=3,
      color="grey",
      help="One point per data to segment",
      showSelected="Split",
      data=set_data_dt)+
    geom_path(aes(
      position, mean, key=1),
      data=seg_path_dt,
      size=1,
      help="Segment means",
      showSelected=c("Split","segments")),
  loss=ggplot()+
    theme_bw()+
    make_tallrect(
      loss_dt, "segments")+
    geom_line(aes(
      segments, Total_Square_Loss, color=set, group=paste(set,Split)),
      showSelected="set",
      clickSelects="Split",
      alpha_off=0.2,
      size=3,
      help="Total squared error for each set",
      data=loss_dt)+
    geom_point(aes(
      segments, Total_Square_Loss,
      key=segments,
      color=set,
      fill=point),
      size=4,
      help="Black point emphasizes min of each curve",
      showSelected=c("Split","set"),
      data=loss_dt)+
    scale_fill_manual(values=c(
      min="black",
      other=NA))+
    scale_color_manual(
      guide="none",
      values=set_colors)+
    scale_y_log10()+
    facet_grid(set ~ ., labeller=label_both, scales="free")+
    scale_x_continuous(),
  duration=list(
    Split=1000,
    segments=1000),
  first=list(segments=5),
  source="https://github.com/tdhock/binseg-model-selection/blob/main/figure-cv-interactive/figure-several-splits.R",
  title="Cross-validation for change-point model selection, several splits"
)
viz
if(FALSE){
  animint2pages(viz, "figure-several-splits", chromote_sleep_seconds=5)
}
