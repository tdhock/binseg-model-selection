source("packages.R")
## are there several optimal trees? N=50,m=3,s=4
node.dt.list <- list()
for(N.data in 60:100){
  min.seg.len <- 5L
  N.segs <- 10
  N.changes <- N.segs-1L
  f.dt <- data.table(d=0:N.changes)[, data.table(
    s=if(N.changes==d)N.data else
      seq(min.seg.len*(d+1), N.data-min.seg.len*(N.changes-d))
  ), by=d]
  setkey(f.dt, d, s)
  g <- function(size){
    ch <- size-min.seg.len*2+1
    ifelse(ch<0, 0, ch)
  }
  tiebreak.list <- list(
    "Equal depth"=function(DT)DT[which.min(d.diff)],
    ##unequal.depth=function(DT)DT[which.max(d.diff)],
    "Equal size"=function(DT)DT[which.min(s.diff)],
    "Unequal size"=function(DT)DT[which.max(s.diff)])
  for(tiebreak in names(tiebreak.list)){
    tiebreak.fun <- tiebreak.list[[tiebreak]]
    for(d.value in 0:N.changes){
      out.dt <- f.dt[J(d.value)]
      out.g <- g(out.dt$s)
      out.f <- if(d.value==0) 0 else {
        cost.dt <- data.table(
          d.under=seq(0, floor((d.value-1)/2))
        )[, {
          d.over <- d.value-d.under-1
          data.table(s.out=out.dt$s)[, {
            seq.end <- min(
              if(d.over == d.under)floor(s.out/2),
              s.out+(d.under-d.value)*min.seg.len)
            seq.start <- (d.under+1)*min.seg.len
            data.table(d.over, s.under=seq.start:seq.end)
          }, by=s.out]
        }, by=d.under]
        cat(sprintf(
          "Ndata=%d Nsegs=%d iteration=%d %d cost values\n", 
          N.data, N.segs, d.value, nrow(cost.dt)))
        cost.dt[, s.over := s.out - s.under]
        cost.dt[, d.diff := abs(d.under-d.over)]
        cost.dt[, s.diff := abs(s.under-s.over)]
        cost.dt[, f.over := f.dt[J(d.over, s.over)]$f]
        cost.dt[, f.under := f.dt[J(d.under, s.under)]$f]
        cost.dt[, f := f.under+f.over]
        cost.dt[, .(s.under,d.under,f,s.over,d.over)]
        best.cost <- cost.dt[out.dt, {
          min.rows <- .SD[f==min(f)]
          tiebreak.fun(min.rows)
        }, keyby=.EACHI, on=.(s.out=s)]
        f.dt[out.dt, `:=`(
          s1=best.cost$s.under, d1=best.cost$d.under,
          s2=best.cost$s.over, d2=best.cost$d.over)]
        best.cost$f
      }
      f.dt[out.dt, f := out.f+out.g]
    }
    tree.dt <- binsegRcpp::get_complexity_best_optimal_tree(f.dt)
    nodes <- binsegRcpp::tree_layout(tree.dt)
    node.dt.list[[paste(N.data, tiebreak)]] <- data.table(
      N.data, tiebreak, nodes)
  }
}
node.dt <- do.call(rbind, node.dt.list)

# hilite for transition slides.
for(N.hilite in 65:75){
  node.hilite <- node.dt[N.data==N.hilite]
  cost <- node.hilite[depth==0]
  gg <- ggplot()+
    ggtitle(sprintf(
      "Min segment length=%d, Splits=%d, Cost=%d",
      min.seg.len, N.changes, cost$f[1]))+
    theme_bw()+
    theme(panel.spacing=grid::unit(0,"lines"))+
    facet_grid(. ~ tiebreak)+
    geom_segment(aes(
      x,depth,xend=parent.x,yend=parent.depth),
      data=node.hilite)+
    geom_label(aes(
      x,depth,label=size),
      size=2.5,
      data=node.hilite)+
    scale_y_reverse("",breaks=NULL)+
    scale_x_continuous("",breaks=NULL)
  png(
    sprintf("figure-optimal-trees-%d.png", N.hilite),
    width=6, height=2, units="in", res=200)
  print(gg)
  dev.off()
}

node.hilite <- node.dt[N.data %in% c(60, 71, 72, 80) & tiebreak=="Equal size"]
gg <- ggplot()+
  theme_bw()+
  theme(panel.spacing=grid::unit(0,"lines"))+
  facet_grid(. ~ N.data, labeller=label_both)+
  geom_segment(aes(
    x,depth,xend=parent.x,yend=parent.depth),
    data=node.hilite)+
  geom_label(aes(
    x,depth,label=size),
    size=2.5,
    data=node.hilite)+
  scale_y_reverse("",breaks=NULL)+
  scale_x_continuous("",breaks=NULL)
png(
  sprintf("figure-optimal-trees-some.png", N.hilite),
  width=8, height=2, units="in", res=200)
print(gg)
dev.off()

gg <- ggplot()+
  ggtitle(sprintf(
    "Min segment length=%d, Splits=%d",
    min.seg.len, N.changes))+
  theme_bw()+
  theme(panel.spacing=grid::unit(0,"lines"))+
  facet_grid(N.data ~ tiebreak)+
  geom_segment(aes(
    x,depth,xend=parent.x,yend=parent.depth),
    data=node.dt)+
  geom_label(aes(
    x,depth,label=size),
    size=2.5,
    data=node.dt)+
  scale_y_reverse("",breaks=NULL)+
  scale_x_continuous("",breaks=NULL)
png(
  "figure-optimal-trees.png", 
  width=6, height=60, units="in", res=200)
print(gg)
dev.off()
