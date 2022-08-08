source("packages.R")
references_funs <- list(
  "N log N"=function(N)log10(N*log(N)),
  "N^2"=function(N)2*log10(N),
  "N^3"=function(N)3*log10(N))
references <- function
(N, empirical, lower.limit=min(empirical),
  fun.list=NULL
){
  fun.latex <- reference <- . <- NULL
  if(is.null(fun.list))fun.list <- references_funs
  data.table(fun.latex=names(fun.list))[, {
    fun <- fun.list[[fun.latex]]
    log10.vec <- fun(N)
    one.fun <- data.table(
      N, empirical,
      reference=10^(log10.vec-max(log10.vec)+max(log10(empirical)))
    )
    above <- one.fun[lower.limit < reference]
    if(1 < nrow(above)){
      above
    }else{
      one.fun[(.N-1):.N]
    }
  }, by=.(fun.latex, fun.name=gsub("\\", "", fun.latex, fixed=TRUE))]
}
references_best <- function(DT, fun.list=NULL){
  unit.col.vec <- c(seconds="seconds_median")
  ref.dt.list <- list()
  metric.dt.list <- list()
  for(unit in names(unit.col.vec)){
    col.name <- unit.col.vec[[unit]]
    values <- DT[[col.name]]
    only.positive <- values[0 < values]
    if(length(only.positive)){
      lower.limit <- min(only.positive)
      all.refs <- DT[
      , references(N.data, .SD[[col.name]], lower.limit, fun.list)
      , by=.(Package, case)]
      all.refs[, rank := rank(N), by=.(Package, case, fun.name)]
      second <- all.refs[rank==2]
      second[, dist := log10(empirical/reference) ]
      second[, sign := sign(dist)]
      l.cols <- list(
        overall=c("Package","case"), 
        each.sign=c("Package","case","sign"))
      for(best.type in names(l.cols)){
        by <- l.cols[[best.type]]
        second[
        , paste0(best.type,".rank") := rank(abs(dist))
        , by=by]
      }
      ref.dt.list[[unit]] <- data.table(unit, all.refs[
        second,
        on=.(Package, case, fun.name, fun.latex)])
      best <- second[overall.rank==1, .(Package, case, fun.name, fun.latex)]
      metric.dt.list[[unit]] <- data.table(unit, DT[
        best, on=.(Package, case)
      ][, `:=`(
        expr.class=paste0(Package, case,"\n",fun.name),
        expr.latex=sprintf("%s %s\n$O(%s)$", Package, case, fun.latex),
        empirical=get(col.name)
      )])
    }
  }
  structure(list(
    references=do.call(rbind, ref.dt.list),
    measurements=do.call(rbind, metric.dt.list)),
    class="references_best")
}
timing.stats <- data.table::fread("figure-timings-data.csv")
timing.stats[, new.pkg := ifelse(
  package %in% names(disp.pkg), disp.pkg[package], package)]
timing.stats[, Package := sub("[.]", "\n", new.pkg, perl=TRUE)]
ref.dt <- rbind(
  data.table(seconds=1, unit="1 second"),
  data.table(seconds=60, unit="1 minute"))
ruptures <- timing.stats[
  package %in% c(
    "ruptures.cumsum","ruptures.l2","changepoint",
    "binsegRcpp.list","binsegRcpp.multiset")]
rbest <- references_best(ruptures)
rank1 <- rbest$references[reference > 1e-3]
gg <- ggplot()+
  theme_bw()+
  theme(
    legend.position="none",
    panel.spacing=grid::unit(0, "lines"))+
  facet_grid(case ~ Package, labeller=label_both)+
  geom_hline(aes(
    yintercept=seconds),
    data=ref.dt,
    color="grey")+
  geom_text(aes(
    10, seconds, label=unit),
    data=ref.dt,
    size=3,
    hjust=0,
    vjust=1.1,
    color="grey50")+
  ## directlabels::geom_dl(aes(
  ##   N.data, seconds_median, color=case, label=case),
  ##   data=ruptures,
  ##   method=list(cex=0.6, "last.polygons"))+
  geom_ribbon(aes(
    N.data, ymin=seconds_min, ymax=seconds_max, fill=case),
    alpha=0.5,
    data=ruptures)+
  geom_line(aes(
    N.data, seconds_median),
    size=1,
    data=ruptures)+
  scale_x_log10(
    "Number of data points to segment (log scale)",
    breaks=10^seq(1, 6, by=1))+
  scale_y_log10(
    "Computation time (seconds, log scale)\nMedian line and min/max band over 5 timings",
    breaks=10^seq(-10,10))+
  geom_line(aes(
    N, reference, group=paste(Package, case, fun.name)),
    color="red",
    data=rank1)+
  directlabels::geom_dl(aes(
    N, reference, 
    label=fun.name, 
    label.group=paste(fun.name)),
    color="red",
    method="bottom.polygons",
    data=rank1)+
  coord_cartesian(ylim=c(1e-4, 1e3))
png("figure-compare-ruptures.png", width=12, height=8, units="in", res=200)
print(gg)
dev.off()
