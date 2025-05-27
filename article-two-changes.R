set.seed(8)
two_changes <- c(rnorm(7, 1), rnorm(10, 3), rnorm(5, 0))
(is.validation.vec <- rep(c(TRUE, FALSE), l=length(two_changes)))
valid_fit <- binsegRcpp::binseg("mean_norm", two_changes, 
  is.validation.vec = is.validation.vec)
valid_fit
(vsegs <- coef(valid_fit, 3L))
vdata <- data.table::data.table(value = two_changes,
 position = seq_along(two_changes), set = ifelse(
   is.validation.vec, "validation", "subtrain"))
valid_join <- vsegs[vdata, .(position, mean, value, set, segments), on=.(start.pos<position, end.pos>position)]
valid_join[, .(loss=sum((mean-value)^2)), by=.(set, segments)]
valid_fit$splits[3, .(segments, loss, validation.loss)]
library(ggplot2)
gg <- ggplot() +
  theme_bw()+
  scale_fill_manual(values = c(subtrain = "black", validation = "white")) +
  scale_color_manual(values = c("validation loss" = "violet", "subtrain mean" = "green")) +
  geom_point(aes(position, value, fill = set), shape = 21, data = vdata) +
  geom_segment(aes(
    start.pos, mean, xend = end.pos, yend = mean, color=line),
    data = data.frame(vsegs, line="subtrain mean")) + 
    scale_x_continuous(breaks=1:22)+
  geom_segment(aes(
    position, mean,
    color=line,
    xend=position, yend=value),
    size=1,
    data=data.frame(valid_join[set=="validation"], line="validation loss"))+
  theme(panel.grid.minor=element_blank())
pdf("article-two-changes.pdf", height=2.1)
print(gg)
dev.off()
