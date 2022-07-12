source("packages.R")
sha.vec <- c(
  "before l1loss"="5a150dfbf4f8934e802b6cead6add3c09c111ace",
  "unordered_map"="bcef04479edcd69377b5aade566f26765da5f6da",
  "mvl_construct"="5942af606641428315b0e63c7da331c4cd44c091",
  "var_param"="2a8dbd5346543bffa2faca293e14735c9ee4e9ee",
  "all l1loss"="71de83db39d14c75d6d85f9d2e133a951d0b0122")
pkg.path <- "~/R/binsegRcpp"
pkg.DESC <- file.path(pkg.path, "DESCRIPTION")
DESC.mat <- read.dcf(pkg.DESC)
Package <- DESC.mat[,"Package"]
new.path <- file.path(tempdir(), Package)
new.DESC <- file.path(new.path, "DESCRIPTION")
unlink(new.path, recursive=TRUE, force=TRUE)
system(paste("cp -r", pkg.path, new.path))
owd <- setwd(new.path)
lib.path <- file.path(owd, "library")
for(commit.name in names(sha.vec)){
  sha <- sha.vec[[commit.name]]
  new.Package <- sprintf("%s.%s", Package, sha)
  system("git checkout -- NAMESPACE DESCRIPTION R/*.R src/*.cpp")
  system(paste("git checkout", sha))
  DESC.mat[,"Package"] <- new.Package
  write.dcf(DESC.mat, new.DESC)
  unlink(file.path(lib.path, new.Package), recursive=TRUE, force=TRUE)
  system(sprintf("sed -i 's/%s/%s/g' R/*.R NAMESPACE", Package, new.Package))
  Rcpp::compileAttributes(new.path)
  system(paste("R CMD INSTALL -l",lib.path,new.path))
}
setwd(owd)
.libPaths("library")
seconds.limit <- 100
do.sub <- function(...){
  mcall <- match.call()
  L <- as.list(mcall[-1])
  for(arg.name in names(L)){
    maybe.lang <- L[[arg.name]]
    if(is.language(maybe.lang)){
      L[[arg.name]] <- substitute(
        result.list[[NAME]] <- EXPR,
        list(NAME=arg.name, EXPR=maybe.lang))
    }
  }
  L
}

timing.dt.list <- list()
done.list <- list()
for(N.data.exp in 2:20){#2^20 = 1,048,576
  N.data <- 2^N.data.exp
  max.segs <- as.integer(N.data/2)
  max.changes <- max.segs-1L
  print(N.data)
  set.seed(1)
  data.list <- list(
    flat=rnorm(N.data),
    best=1:N.data,
    worst=rep(c(0,1),l=N.data))
  for(case in names(data.list)){
    data.vec <- data.list[[case]]
    cum.data.vec <- cumsum(c(0,data.vec))
    data.mat <- matrix(data.vec)
    result.list <- list()
    m.args <- list(times=5)
    for(commit.name in names(sha.vec)){
      sha <- sha.vec[[commit.name]]
      new.Package <- sprintf("%s.%s", Package, sha)
      binseg.symbol <- str2lang(paste0(new.Package, "::binseg"))
      m.args[[commit.name]] <- substitute(BINSEG(
        "mean_norm",data.vec, max.segs, container.str="multiset"),
        list(BINSEG=binseg.symbol))
    }
    m.args[names(done.list[[case]])] <- NULL
    if(length(m.args) > 1){
      N.df <- do.call(microbenchmark::microbenchmark, m.args)
      N.dt <- data.table(N.df)
      N.dt[, seconds := time/1e9]
      N.dt[, package := paste(expr)]
      N.stats <- data.table::dcast(
        N.dt,
        package ~ .,
        list(median, min, max, timings=length),
        value.var="seconds")
      done.pkgs <- N.stats[seconds_median > seconds.limit, paste(package)]
      done.list[[case]][done.pkgs] <- TRUE
      timing.dt.list[[paste(N.data, case)]] <- data.table(
        N.data, max.segs, case, N.stats)
    }
  }
}
timing.dt <- do.call(rbind, timing.dt.list)

data.table::fwrite(timing.dt, "figure-timings-versions-data.csv")
