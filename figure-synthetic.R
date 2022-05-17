a <- 1
x <- sqrt(24/9)
e <- 2
b <- 10
linear <- c(b+e+x,b+e,b-e,b-e-x)
m2 <- c(b+e+x/2, b+e+x/2, b-e-x/2, b-e-x/2)
m1 <- mean(linear)
loss <- function(m)sum((linear-m)^2)
loss(m1)-loss(m2)
rbind(loss(m1), 2*e^2+2*(e+x)^2)
rbind(loss(m2), x^2)
data.vec <- c(a,-a,a,-a,linear)
plot(data.vec)
fit <- binsegRcpp::binseg_normal(data.vec)
library(data.table)
some.sizes <- 2:6
decrease.dt <- fit$splits[, data.table(
  loss.decrease=diff(loss),
  end=end[-1],
  segments=seq(2, length(data.vec)))
  ][J(some.sizes), on="segments"]
some.segs <- coef(fit, some.sizes)
synth.dt <- data.table(data.value=data.vec, position=seq_along(data.vec))
library(ggplot2)
model.color <- "red"
decrease.dt[, split := segments-1]
some.segs[, split := segments-1]
gg <- ggplot()+
  facet_grid(. ~ split, labeller=label_both)+
  geom_point(aes(
    position, data.value),
    size=2,
    data=synth.dt)+
  geom_segment(aes(
    start-0.5, mean,
    xend=end+0.5, yend=mean),
    color=model.color,
    size=1,
    data=some.segs)+
  geom_vline(aes(
    xintercept=end+0.5),
    color=model.color,
    data=decrease.dt)+
  geom_text(aes(
    end+0.5, Inf,
    label=sprintf("loss decrease\n=%.2f", -loss.decrease),
    hjust=ifelse(end<5, 0, 1)),
    color="red",
    vjust=1,
    size=2.5,
    data=decrease.dt)+
  scale_y_continuous(limits=c(NA, 18))+
  scale_x_continuous(breaks=seq(1, length(data.vec)))+
  theme_bw()+
  theme(panel.spacing=grid::unit(0,'lines'))
png("figure-synthetic.png", width=7, height=2, units="in", res=200)
print(gg)
dev.off()
