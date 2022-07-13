source("packages.R")
sha.vec <- c(
  cv="908b77c411bc7f4fcbcf53759245e738ae724c3e",
  ##"unordered_map"="bcef04479edcd69377b5aade566f26765da5f6da",
  ##"mvl_construct"="5942af606641428315b0e63c7da331c4cd44c091",
  ##"new dist meth"="544a4df46ba53ccae8a82a3c67911007cd008000",
  ##"laplace tests"="61ad58075ae0399d69a8565949ff7aa944f0a461",
  ## "vloss tests"="53ab142f5b866e80a2d8581f694499f3531c47c1",
  ## "std::vector"="4732a3af89776461f0398892221c1d09236cda1f",
  ## "consolidate"="03af3dcec278118808bad5f1bf4dc803f28030e8",
  ## "simplify con"="2352f6efc0c11fac94e802c960b9ea38a30e0539",
  ## "var_param"="2a8dbd5346543bffa2faca293e14735c9ee4e9ee",
  ## "all l1loss"="71de83db39d14c75d6d85f9d2e133a951d0b0122",
  ## "candidate_ptr arg"="94277b39558fc8c51da2689a738e1b990623bb68",
  ## "Split pointer arg"="0200e7c55603707b2393bcaf29f1ca64250a3986",
  ## "ParamsLoss ptr"="b77a6b8f3b74da4e0f3eb1c0b99cf9dfa310ab11",
  ##"mean var ptr"="26279384340413007c0839b8ba8652b0f6c59f1e",
  "rm unord map"="dcd0808f52b0b9858352106cc7852e36d7f5b15d",
  "no var est"="28f2e5aa56a5c7e4147fc65a86a4f4936cecaa03")
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
dir.create(lib.path, showWarnings = FALSE)
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
seconds.limit <- 0.1
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
      binseg.symbol <- str2lang(paste0(new.Package, "::binseg_normal"))
      m.args[[commit.name]] <- substitute(BINSEG(
        data.vec,max.segs),
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
