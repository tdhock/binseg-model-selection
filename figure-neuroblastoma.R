source("packages.R")
Sys.setenv(RETICULATE_PYTHON=if(.Platform$OS.type=="unix")
  "/home/tdhock/.local/share/r-miniconda/envs/cs570s22/bin/python"
  else "~/Miniconda3/envs/cs570s22/python.exe")
reticulate::use_condaenv("cs570s22", required=TRUE)
ruptures <- reticulate::import("ruptures")
data(neuroblastoma,package="neuroblastoma")
nb.profiles <- data.table(neuroblastoma[["profiles"]])
one.pid.chr <- nb.profiles[profile.id==2 & chromosome==2]
one.pid.chr[, data.i := .I]
data.vec <- one.pid.chr$logratio
gg <- ggplot()+
  scale_x_continuous(
    "Position/index in data sequence",
    limits=c(0, nrow(one.pid.chr)+1))+
  scale_y_continuous(
    "DNA copy number logratio to segment")+
  geom_point(aes(
    data.i, logratio),
    shape=21,
    data=one.pid.chr)
gg

N.data <- length(data.vec)
cum.data.vec <- cumsum(c(0,data.vec))
data.mat <- matrix(data.vec)
max.changes <- 4L
max.segs <- max.changes+1L
binseg_instance <- ruptures$Binseg(min_size=1L, jump=1L)
cpt.fit <- changepoint::cpt.mean(data.vec, method="BinSeg", Q=max.changes)
wbs.fit <- wbs::sbs(data.vec)
wbs.dt <- data.table(wbs.fit$res)[order(-min.th, scale)]
binseg.fit <- binsegRcpp::binseg_normal(data.vec, max.segs)
loss.dt.list <- list()
seg.dt.list <- list()
last.ruptures <- N.data
ruptures.ord <- list()
for(n.changes in 1:max.changes){
  n.segs <- n.changes+1L
  end.list <- list(
    ruptures=binseg_instance$fit(data.mat)$predict(n_bkps=n.changes),
    changepoint=c(N.data,cpt.fit@cpts.full[n.changes,]),
    ##wbs=wbs.dt[, c(N.data, cpt)][1:n.segs],
    "binsegRcpp=wbs"=coef(binseg.fit, n.segs)$end
  )
  ruptures.ord[[n.changes]] <- with(
    end.list, ruptures[!ruptures %in% last.ruptures])
  last.ruptures <- end.list$ruptures
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

loss.wide <- dcast(
  loss.dt,
  package ~ n.changes,
  value.var="total.square.loss"
)[c("binsegRcpp=wbs","ruptures","changepoint"), on="package"]
loss.wide.xt <- xtable(loss.wide)
print(loss.wide.xt, floating=FALSE, include.rownames=FALSE)

seg.dt <- do.call(rbind, seg.dt.list)
split.dt.list <- list()
change.list <- list(
  changepoint=end.list$changepoint[-1],
  ruptures=as.integer(ruptures.ord),
  ##wbs=end.list$wbs[-1],
  "binsegRcpp=wbs"=binseg.fit$splits$end[-1])
package.y <- seq(2, 4, l=length(change.list))
names(package.y) <- names(change.list)
for(package in names(change.list)){
  end <- change.list[[package]]
  split.dt.list[[package]] <- data.table(
    package, package.y=package.y[[package]],
    split.i=seq_along(end), end)
}
split.dt <- do.call(rbind, split.dt.list)
model.color <- "red"
gg+
  geom_segment(aes(
    start-0.5, mean,
    xend=end+0.5, yend=mean),
    data=seg.dt,
    color=model.color)+
  geom_vline(aes(
    xintercept=start-0.5),
    data=seg.dt[1 < start],
    color=model.color)+
  facet_grid(package ~ n.changes)
pkg.dt <- data.table(package=names(package.y), package.y)
pkg.loss <- loss.dt[, .(
  loss.values=paste(sprintf("%.2f", total.square.loss), collapse=", ")
), by=package][pkg.dt, on="package"]
out <- gg+
  theme_bw()+
  geom_vline(aes(
    xintercept=end+0.5),
    data=split.dt,
    color="grey50")+
  geom_text(aes(
    65, package.y, label=paste(package, "loss values=", loss.values)),
    hjust=1,
    data=pkg.loss)+
  geom_label(aes(
    end+0.5, package.y, label=split.i),
    label.padding = unit(0.1, "lines"),
    size=3,
    data=split.dt)+
  coord_cartesian(
    xlim=c(15,75))
png("figure-neuroblastoma.png", width=7, height=3, units="in", res=200)
print(out)
dev.off()
