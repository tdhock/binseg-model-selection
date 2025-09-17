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

selection_dt <- rbindlist(selection_dt_list)
selection_cum <- data.table(num_splits=seq(4, total.validation.sets, by=n.folds))[, {
  selection_dt[validation_set_i<=num_splits, .(count=.N), by=.(split_type,segments)]
}, by=num_splits][
, percent := 100*count/sum(count), by=.(num_splits,split_type)][
, tooltip := sprintf(
  "%d segments selected in %d/%d=%.1f%% validation sets, split type %s",
  segments, count, num_splits, percent, split_type)
]
selection_cum[num_splits==total.validation.sets][order(-count)]
selection_most_freq <- selection_cum[, .SD[count==max(count)], by=.(num_splits,split_type)]

animint(
  ggplot()+
    theme_bw()+
    geom_tile(aes(
      num_splits, segments,
      tooltip=tooltip,
      fill=percent),
      color=NA,
      data=selection_cum)+
    scale_color_manual(
      values=c(largest="black"))+
    geom_point(aes(
      num_splits, segments,
      tooltip=tooltip,
      color=frequency),
      fill="transparent",
      data=data.table(frequency="largest", selection_most_freq))+
    scale_fill_gradient(low="white", high="red")+
  facet_grid(split_type ~ ., labeller=label_both)
)

loss_dt <- rbindlist(loss_dt_list)[, Total_Square_Loss := ifelse(loss<1e-10, 0, loss)][]
valid_dt <- loss_dt[set=="validation"][
, log10.loss := log10(loss)][
, norm.valid.loss := (log10.loss-min(log10.loss))/(max(log10.loss)-min(log10.loss)), by=.(split_type, validation_set_i)]
valid_min <- valid_dt[point=="min"]
valid_min_stats <- dcast(
  valid_min,
  split_type + validation_set_i ~ .,
  list(min, max),
  value.var="segments"
)[order(segments_min, segments_max)][
, split_sorted_index := 1:.N, by=split_type]
i_to_sorted <- valid_min_stats[, .(split_type, split_sorted_index, validation_set_i)]
add_sorted <- function(DT)i_to_sorted[DT, on=.(split_type, validation_set_i)][
, Split := paste(split_type, split_sorted_index)]
valid_tiles <- add_sorted(valid_dt)
tallrect_dt <- valid_tiles[segments==1]
my_rbind <- function(L)add_sorted(rbindlist(L))
(set_data_dt <- my_rbind(set_data_dt_list))
library(animint2)
dput(RColorBrewer::brewer.pal(Inf, "Set2"))
c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", 
  "#E5C494", "#B3B3B3")
set_colors <- c(
  subtrain="#66C2A5",
  validation="#FC8D62")
(residual_dt <- my_rbind(residual_dt_list))
seg_path_dt <- data.table(offset=0.5*c(1,-1))[, residual_dt[, .(
  position=position-offset, mean, segments, Split, data_i=position
)], by=offset][order(segments, data_i, position)]
loss_sorted <- add_sorted(loss_dt)
viz <- animint(
  splits=ggplot()+
    ggtitle("Select subtrain/validation split")+
    theme_bw()+
    theme_animint(
      width=700, colspan=2, last_in_row=TRUE)+
    geom_tile(aes(
      split_sorted_index, segments,
      fill=norm.valid.loss),
      data=valid_tiles,
      color=NA)+
    geom_point(aes(
      split_sorted_index, segments,
      color=point),
      fill="transparent",
      data=add_sorted(valid_min))+
    geom_tallrect(aes(
      xmin=split_sorted_index-0.5,
      xmax=split_sorted_index+0.5,
      ymin=-Inf, ymax=Inf),
      data=tallrect_dt,
      clickSelects="Split",
      fill="blue",
      alpha=0.5,
      color=NA)+
    scale_color_manual(values=c(min="violet"))+
    scale_fill_gradient(low="white", high="black")+
    facet_grid(split_type ~ ., labeller=label_both),
  loss=ggplot()+
    ggtitle("Select segments")+
    theme_bw()+
    theme_animint(width=300, height=300)+
    geom_line(aes(
      segments, Total_Square_Loss, color=set, group=set,
      key=set),
      showSelected=c("set","Split"),
      size=3,
      help="Total squared error for each set",
      data=loss_sorted)+
    geom_point(aes(
      segments, Total_Square_Loss,
      key=segments,
      color=set,
      fill=point),
      size=4,
      help="Black point emphasizes min of each curve",
      showSelected=c("Split","set"),
      data=loss_sorted)+
    make_tallrect(
      loss_dt, "segments")+
    scale_fill_manual(values=c(
      min="black",
      other=NA))+
    scale_color_manual(
      values=set_colors)+
    scale_y_log10()+
    facet_grid(set ~ ., labeller=label_both, scales="free")+
    scale_x_continuous(),
  data=ggplot()+
    ggtitle("Residual errors of selected model")+
    theme_bw()+
    theme_animint(
      width=400, height=300)+
    theme(legend.position="none")+
    scale_fill_manual(values=set_colors)+
    scale_color_manual(values=set_colors)+
    geom_segment(aes(
      position, signal, key=position,
      color=set,
      xend=position, yend=mean),
      showSelected=c("set","Split","segments"),
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
      showSelected=c("set","Split"),
      data=set_data_dt)+
    geom_path(aes(
      position, mean, key=1),
      data=seg_path_dt,
      size=1,
      help="Segment means",
      showSelected=c("Split","segments")),
  duration=list(
    Split=1000,
    segments=1000),
  first=list(segments=5),
  source="https://github.com/tdhock/binseg-model-selection/blob/main/figure-cv-interactive/figure-most-frequently-selected.R",
  title="Most frequently selected change-point model using cross-validation")

viz
if(FALSE){
  animint2pages(viz, "figure-binseg-cv-most-frequently-selected", chromote_sleep_seconds=5)
}
