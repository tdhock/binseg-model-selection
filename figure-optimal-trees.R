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
    ##decoding.
    level <- 0
    new.id <- 1
    new.nodes <- data.table(
      f.dt[.N, .(d,s,parent.x=NA,parent.y=NA)])
    while(nrow(new.nodes)){
      node.info <- f.dt[new.nodes][order(parent.x,s)]
      node.info[, id := seq(new.id, new.id+.N-1)]
      node.info[, ord := 1:.N]
      node.info[, y := -level]
      node.info[, x := ord-(.N+1)/2]
      new.id <- node.info[.N, new.id+1]
      node.dt.list[[paste(N.data, N.segs, tiebreak, level)]] <- 
        node.info[, data.table(
          N.data, N.segs, tiebreak, level, f,
          ord,id,d,s,x,y,parent.x,parent.y)]
      level <- level+1
      new.nodes <- node.info[, data.table(
        d=c(d1,d2),s=c(s1,s2),parent.x=x,parent.y=y
      )][!is.na(d)][order(parent.x)]
    }
  }
}
node.dt <- do.call(rbind, node.dt.list)

# hilite for transition slides.
for(N.hilite in 65:75){
  node.hilite <- node.dt[N.data==N.hilite]
  cost <- node.hilite[level==0]
  gg <- ggplot()+
    ggtitle(sprintf(
      "Min segment length=%d, Splits=%d, Cost=%d",
      min.seg.len, N.changes, cost$f[1]))+
    theme_bw()+
    theme(panel.spacing=grid::unit(0,"lines"))+
    facet_grid(. ~ tiebreak)+
    geom_segment(aes(
      x,y,xend=parent.x,yend=parent.y),
      data=node.hilite)+
    geom_label(aes(
      x,y,label=s),
      size=3,
      data=node.hilite)+
    scale_y_continuous("",breaks=NULL)+
    scale_x_continuous("",breaks=NULL)
  png(
    sprintf("figure-optimal-trees-%d.png", N.hilite),
    width=6, height=3, units="in", res=200)
  print(gg)
  dev.off()
}

gg <- ggplot()+
  ggtitle(sprintf(
    "Min segment length=%d, Splits=%d",
    min.seg.len, N.changes))+
  theme_bw()+
  theme(panel.spacing=grid::unit(0,"lines"))+
  facet_grid(N.data ~ tiebreak)+
  geom_segment(aes(
    x,y,xend=parent.x,yend=parent.y),
    data=node.dt)+
  geom_label(aes(
    x,y,label=s),
    size=3,
    data=node.dt)+
  scale_y_continuous("",breaks=NULL)+
  scale_x_continuous("",breaks=NULL)
png(
  "figure-optimal-trees-all.png", 
  width=6, height=60, units="in", res=200)
print(gg)
dev.off()
