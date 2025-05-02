changepoint::cpt.mean(1:2, method="BinSeg", Q=2)@cpts.full
changepoint::cpt.mean(1:3, method="BinSeg", Q=2)@cpts.full
changepoint::cpt.mean(1:4, method="BinSeg", Q=3)@cpts.full
changepoint::cpt.mean(1:5, method="BinSeg", Q=3)@cpts.full
changepoint::cpt.mean(1:6, method="BinSeg", Q=4)@cpts.full
changepoint::cpt.mean(1:8, method="BinSeg", Q=3)@cpts.full

library(data.table)
data(neuroblastoma,package="neuroblastoma")
nb.profiles <- data.table(neuroblastoma[["profiles"]])
one.pid.chr <- nb.profiles[profile.id==2 & chromosome==2]
one.pid.chr[, data.i := .I]
data.vec <- one.pid.chr$logratio
N.data <- length(data.vec)
cum.data.vec <- cumsum(c(0,data.vec))
max.segs <- 5
max.changes <- max.segs-1
binseg.fit <- binsegRcpp::binseg_normal(data.vec, max.segs)
cpt.fit <- changepoint::cpt.mean(data.vec, method="BinSeg", Q=max.changes)
loss.dt.list <- list()
seg.dt.list <- list()
for(n.changes in 1:max.changes){
  n.segs <- n.changes+1L
  end.list <- list(
    changepoint=c(N.data, cpt.fit@cpts.full[n.changes,]),
    binsegRcpp=coef(binseg.fit, n.segs)$end
  )
  for(package in names(end.list)){
    end <- sort(end.list[[package]])
    seg.size <- diff(c(0,end))
    seg.mean <- (cum.data.vec[end+1]-cum.data.vec[end-seg.size+1])/seg.size
    data.mean <- rep(seg.mean, seg.size)
    seg.dt.list[[paste(n.changes, package)]] <- data.table(
      n.changes,
      package,
      start=end-seg.size+1L,
      end,
      mean=seg.mean)
    loss.dt.list[[paste(n.changes, package)]] <- data.table(
      n.changes, package, total.square.loss=sum((data.mean-data.vec)^2))
  }
}
loss.dt <- do.call(rbind, loss.dt.list)
dcast(loss.dt, package ~ n.changes, value.var="total.square.loss")
seg.dt <- do.call(rbind, seg.dt.list)
loss.y <- c(
  binsegRcpp=6,
  changepoint=4)
library(ggplot2)
gg <- ggplot()+
  theme_bw()+
  scale_color_manual(values=c(binsegRcpp="red", changepoint="black"))+
  scale_size_manual(values=c(binsegRcpp=1, changepoint=2))+
  geom_text(aes(
    40, loss.y[package], color=package,
    label=sprintf(
      "%s (%s) loss= %.2f",
      package,
      ifelse(package=="binsegRcpp", "correct", "incorrect"),
      total.square.loss)),
    hjust=1,
    data=loss.dt)+
  scale_x_continuous(
    "Position/index in data sequence",
    breaks=seq(0,1000,by=2),
    limits=c(0, nrow(one.pid.chr)+1))+
  scale_y_continuous(
    "DNA copy number logratio to segment")+
  geom_segment(aes(
    start-0.5, mean,
    size=package,
    color=package,
    xend=end+0.5, yend=mean),
    data=seg.dt)+
  geom_vline(aes(
    xintercept=start-0.5,
    size=package,
    color=package),
    linetype="dashed",
    data=seg.dt[1 < start])+
  facet_grid(n.changes ~ ., labeller=label_both)+
  geom_point(aes(
    data.i, logratio),
    color="grey50",
    shape=21,
    fill="white",
    data=one.pid.chr)+
  coord_cartesian(
    expand=FALSE,
    ylim=c(-1, 7),
    xlim=c(19,71))
png("changepoint.bug.png", width=10, height=5, units="in", res=200)
print(gg)
dev.off()
