library(data.table)
library(ggplot2)
pkg.colors <- c(
  "changepoint\narray"="#E41A1C",
  "ruptures\nLRU cache"="#377EB8",
  "ruptures\nl2"="#377EB8",
  "ruptures\ncumsum"="#377EB8",
  "blockcpd\nheap"="#984EA3", 
  "fpop::multiBinSeg\nheap"="#FF7F00", 
  "wbs::sbs\nrecursion"="#A65628",
  "binsegRcpp\nmultiset"="black",
  "binsegRcpp\nlist"="grey50")
disp.pkg <- c(
  blockcpd="blockcpd.heap",
  changepoint="changepoint.array",
  ruptures="ruptures.LRU cache",
  "wbs::sbs"="wbs::sbs.recursion",
  "fpop::multiBinSeg"="fpop::multiBinSeg.heap")

