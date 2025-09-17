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
    Méthode_Division <- sprintf("VC à %d blocs", n.folds)
    for(cv_split_seed in 1:n.seeds){
      set.seed(cv_split_seed)
      fold.vec <- sample(rep(uniq.folds, l=n.data))
      for(valid.fold in uniq.folds){
        is.valid.vec.list[[Méthode_Division]][[paste(cv_split_seed, valid.fold)]] <-
          fold.vec==valid.fold
      }
    }
  }
  for(prop.valid in prop.valid.vec){
    Méthode_Division <- sprintf("%d%% validation", 100*prop.valid)
    prop.vec <- c("sous-entraînement"=1-prop.valid, validation=prop.valid)
    for(split_seed in 1:total.validation.sets){
      set.seed(split_seed)
      is.valid.vec.list[[Méthode_Division]][[split_seed]] <- binsegRcpp::random_set_vec(
        n.data, prop.vec) == "validation"
    }
  }
  for(Méthode_Division in names(is.valid.vec.list)){
    for(validation_set_i in seq_len(total.validation.sets)){
      is.valid <- is.valid.vec.list[[Méthode_Division]][[validation_set_i]]
      set(split_data, j="ensemble", value=ifelse(is.valid, "validation", "sous-entraînement"))
      set_data_dt_list[[paste(
        Méthode_Division, validation_set_i
      )]] <- data.table(
        Méthode_Division, validation_set_i,
        split_data)
      bs.model <- split_data[, binsegRcpp::binseg_normal(
        signal, is.validation.vec=is.valid)]
      coef_dt <- coef(bs.model, 1:nrow(bs.model$splits))
      setkey(coef_dt, start.pos, end.pos)
      over_dt <- foverlaps(coef_dt, split_data)
      residual_dt_list[[paste(
        data_seed, Méthode_Division, validation_set_i
      )]] <- over_dt[, .(
        data_seed, Méthode_Division, validation_set_i,
        position, segments, signal, mean, ensemble
      )]
      loss_dt_list[[paste(
        data_seed, Méthode_Division, validation_set_i
      )]] <- data.table(
        data_seed, Méthode_Division, validation_set_i,
        bs.model$splits[, rbind(
          data.table(ensemble="sous-entraînement", segments, loss),
          data.table(ensemble="validation", segments, loss=validation.loss)
        )]
      )[, point := ifelse(loss==min(loss), "minimum", "autre"), by=ensemble][]
      selection_dt_list[[paste(
        data_seed, Méthode_Division, validation_set_i
      )]] <- data.table(
        data_seed, Méthode_Division, validation_set_i,
        segments=bs.model$splits[which.min(validation.loss), segments]
      )
    }
  }
}

selection_dt <- rbindlist(selection_dt_list)
selection_cum <- data.table(num_splits=seq(4, total.validation.sets, by=n.folds))[, {
  selection_dt[validation_set_i<=num_splits, .(count=.N), by=.(Méthode_Division,segments)]
}, by=num_splits][
, percent := 100*count/sum(count), by=.(num_splits,Méthode_Division)][
, tooltip := sprintf(
  "%d segments sélectionnés dans %d/%d=%.1f%% ensembles validation, méthode division : %s",
  segments, count, num_splits, percent, Méthode_Division)
]
selection_cum[num_splits==total.validation.sets][order(-count)]
selection_most_freq <- selection_cum[, .SD[count==max(count)], by=.(num_splits,Méthode_Division)]
animint(
  ggplot()+
    theme_bw()+
    theme_animint(width=600, height=600)+
    geom_tile(aes(
      segments, num_splits,
      tooltip=tooltip,
      fill=percent),
      color=NA,
      data=selection_cum)+
    scale_color_manual(
      values=c(largest="black"))+
    geom_point(aes(
      segments, num_splits,
      tooltip=tooltip,
      color=frequency),
      fill="transparent",
      data=data.table(frequency="largest", selection_most_freq))+
    scale_fill_gradient(low="white", high="red")+
  facet_grid(. ~ Méthode_Division, labeller=label_both)
)

loss_dt <- rbindlist(loss_dt_list)[, Erreur_Quadratique := ifelse(loss<1e-10, 0, loss)][]
valid_dt <- loss_dt[ensemble=="validation"][
, log10.loss := log10(loss)][
, Err_valid_norm := (log10.loss-min(log10.loss))/(max(log10.loss)-min(log10.loss)), by=.(Méthode_Division, validation_set_i)]
valid_min <- valid_dt[point=="minimum"]
valid_min_stats <- dcast(
  valid_min,
  Méthode_Division + validation_set_i ~ .,
  list(min, max),
  value.var="segments"
)[order(segments_min, segments_max)][
, numéro_division_trié := 1:.N, by=Méthode_Division]
i_to_sorted <- valid_min_stats[, .(Méthode_Division, numéro_division_trié, validation_set_i)]
add_sorted <- function(DT)i_to_sorted[DT, on=.(Méthode_Division, validation_set_i)][
, Division := paste(Méthode_Division, numéro_division_trié)]
valid_tiles <- add_sorted(valid_dt)
tallrect_dt <- valid_tiles[segments==1]
my_rbind <- function(L)add_sorted(rbindlist(L))
(set_data_dt <- my_rbind(set_data_dt_list))
dput(RColorBrewer::brewer.pal(Inf, "Set2"))
c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", 
  "#E5C494", "#B3B3B3")
ensemble_colors <- c(
  "sous-entraînement"="#66C2A5",
  validation="#FC8D62")
(residual_dt <- my_rbind(residual_dt_list))
seg_path_dt <- data.table(offset=0.5*c(1,-1))[, residual_dt[, .(
  position=position-offset, mean, segments, Division, data_i=position
)], by=offset][order(segments, data_i, position)]
loss_sorted <- add_sorted(loss_dt)
viz <- animint(
  splits=ggplot()+
    ggtitle("Choisir division sous-entraînement / validation")+
    theme_bw()+
    theme_animint(
      width=600, height=900, rowspan=2)+
    geom_tile(aes(
      segments, numéro_division_trié,
      fill=Err_valid_norm),
      data=valid_tiles,
      color=NA)+
    geom_point(aes(
      segments, numéro_division_trié,
      color=point),
      fill="transparent",
      data=add_sorted(valid_min))+
    geom_widerect(aes(
      ymin=numéro_division_trié-0.5,
      ymax=numéro_division_trié+0.5),
      data=tallrect_dt,
      clickSelects="Division",
      fill="blue",
      alpha=0.5,
      color=NA)+
    scale_color_manual(values=c(minimum="violet"))+
    scale_fill_gradient(low="white", high="black")+
    facet_grid(. ~ Méthode_Division, labeller=label_both),
  loss=ggplot()+
    ggtitle("Choisir segments")+
    theme_bw()+
    theme_animint(width=300, height=550, last_in_row=TRUE)+
    geom_line(aes(
      segments, Erreur_Quadratique, color=ensemble, group=ensemble,
      key=ensemble),
      showSelected=c("ensemble","Division"),
      size=3,
      help="Erreur quadratique pour chaque ensemble",
      data=loss_sorted)+
    geom_point(aes(
      segments, Erreur_Quadratique,
      key=segments,
      color=ensemble,
      fill=point),
      size=4,
      help="Points noirs pour les minima de chaque courbe",
      showSelected=c("Division","ensemble"),
      data=loss_sorted)+
    make_tallrect(
      loss_dt, "segments")+
    scale_fill_manual(values=c(
      minimum="black",
      autre=NA))+
    scale_color_manual(
      values=ensemble_colors)+
    scale_y_log10()+
    facet_grid(ensemble ~ ., labeller=label_both, scales="free")+
    scale_x_continuous(),
  data=ggplot()+
    ggtitle("Erreurs pour le modèle choisi")+
    theme_bw()+
    theme_animint(
      width=400, height=300)+
    theme(legend.position="none")+
    scale_fill_manual(values=ensemble_colors)+
    scale_color_manual(values=ensemble_colors)+
    geom_segment(aes(
      position, signal, key=position,
      color=ensemble,
      xend=position, yend=mean),
      showSelected=c("ensemble","Division","segments"),
      size=1,
      help="Segments verticales pour la différence entre signal et moyenne",
      data=residual_dt)+
    geom_point(aes(
      position, signal,
      fill=ensemble,
      key=position),
      size=3,
      color="grey",
      help="One point per data to segment",
      showSelected=c("ensemble","Division"),
      data=set_data_dt)+
    geom_path(aes(
      position, mean, key=1),
      data=seg_path_dt,
      size=1,
      help="Segment means",
      showSelected=c("Division","segments")),
  duration=list(
    Division=1000,
    segments=1000),
  first=list(segments=5),
  source="https://github.com/tdhock/binseg-model-selection/blob/main/figure-cv-interactive/figure-most-frequently-selected.R",
  out.dir="figure-most-frequently-selected-fr",
  title="Most frequently selected change-point model using cross-validation")

viz
if(FALSE){
  animint2pages(viz, "figure-binseg-cv-most-frequently-selected-fr", chromote_sleep_seconds=5)
}
