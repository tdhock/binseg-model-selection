sim.result <- data.table::fread("figure-select-segments-data.csv")
sim.err <- sim.result[, .(
  zero.one.loss=sum(selected != n.segments),
  L1.loss=sum(abs(selected-n.segments)),
  L2.loss=sum((selected-n.segments)^2)
), by=.(n.segments, method)]
method.dt <- sim.err[, .(
  total=sum(zero.one.loss)
), by=method][order(-total)]
sim.err[, Method := factor(method, method.dt$method)]

library(ggplot2)
gg <- ggplot()+
  geom_tile(aes(
    factor(n.segments), Method,
    fill=zero.one.loss),
    data=sim.err)+
  geom_text(aes(
    factor(n.segments), Method,
    label=zero.one.loss),
    data=sim.err)+
  scale_fill_gradient(
    "Simulations
with
incorrectly
selected
model
size",
low="white", high="red")+
  scale_x_discrete(
    "True number of segments in simulation")+
  scale_y_discrete(
    "Model selection method")+
  theme_bw()
png("figure-select-segments.png", width=6, height=6, units="in", res=200)
print(gg)
dev.off()
