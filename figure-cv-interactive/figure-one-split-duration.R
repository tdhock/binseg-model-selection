library(data.table)
n.segments <- 5
set.seed(1)
(seg.mean.vec <- rnorm(n.segments))
data.per.segment <- 10
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
one_split <- function(L)rbindlist(L)[split_type=="2-fold" & validation_set_i==1]
one_split(selection_dt_list)

loss_dt <- one_split(loss_dt_list)[, Total_Square_Loss := ifelse(loss<1e-10, 0, loss)][]
loss_lab <- loss_dt[segments==max(segments)]
(set_data_dt <- one_split(set_data_dt_list))
library(animint2)
dput(RColorBrewer::brewer.pal(Inf, "Set2"))
c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", 
  "#E5C494", "#B3B3B3")
set_colors <- c(
  subtrain="#66C2A5",
  validation="#FC8D62")
(residual_dt <- one_split(residual_dt_list))
## double check residuals agree with loss.
data.table(
  loss_dt[, .(segments,set,loss)],
  residual_dt[, .(
    RSS=sum((mean-signal)^2)
  ), keyby=.(set, segments)][, .(RSS)]
)[, let(diff=loss-RSS)][]
seg_path_dt <- rbind(
  residual_dt[, .(position=position-0.5, mean, segments, data_i=position)],
  residual_dt[, .(position=position+0.5, mean, segments, data_i=position)]
)[order(segments, data_i, position)]
viz <- animint(
  data=ggplot()+
    theme_bw()+
    scale_fill_manual(values=set_colors)+
    scale_color_manual(values=set_colors)+
    geom_segment(aes(
      position, signal, key=position,
      color=set,
      xend=position, yend=mean),
      showSelected="segments",
      size=1,
      help="Vertical segments show residuals (difference between data and mean)",
      data=residual_dt)+
    geom_point(aes(
      position, signal, fill=set),
      size=3,
      color="grey",
      help="One point per data to segment",
      data=set_data_dt)+
    geom_path(aes(
      position, mean, key=1),
      data=seg_path_dt,
      size=1,
      help="Segment means",
      showSelected="segments"),
  loss=ggplot()+
    theme_bw()+
    geom_line(aes(
      segments, Total_Square_Loss, color=set, group=set),
      showSelected="set",
      help="Total squared error for each set",
      data=loss_dt)+
    geom_point(aes(
      segments, Total_Square_Loss, color=set, fill=point),
      size=4,
      help="Black point emphasizes min of each curve",
      showSelected="set",
      data=loss_dt)+
    scale_fill_manual(values=c(
      min="black",
      other=NA))+
    scale_color_manual(
      guide="none",
      values=set_colors)+
    geom_label_aligned(aes(
      segments, Total_Square_Loss, color=set, label=set),
      hjust=0,
      clickSelects="set",
      data=loss_lab)+
    make_tallrect(
      loss_dt, "segments")+
    scale_y_log10()+
    scale_x_continuous(
      breaks=c(1, seq(5, 25, by=5)),
      limits=c(0, 30)),
  duration=list(segments=1000),
  first=list(segments=5),
  source="https://github.com/tdhock/binseg-model-selection/blob/main/figure-cv-interactive/figure-one-split-duration.R",
  title="Cross-validation for change-point model selection"
)
viz
if(FALSE){
  animint2pages(viz, "figure-binseg-cv-one-split-duration", chromote_sleep_seconds=5)
}
