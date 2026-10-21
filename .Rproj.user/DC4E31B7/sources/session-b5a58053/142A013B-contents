omphc <- function(x, skel = NULL, R = NULL, method = "pvalue", tol = 0.05, algo = "hc",
                   tabu = 10, restart = 10, score = "bic-g", blacklist = NULL, whitelist = NULL) {
  ## score -> "loglik-g", "aic-g", "bic-g", "bge"
  runtime <- proc.time()
  if ( is.null(skel) )  skel <- omphc::omp.network(x, R = R, method = method, tol = tol)
  nama <- colnames(skel$G)
  if ( is.null(nama) )  nama <- paste("X", 1:dim(x)[2], sep = "")
  colnames(x) <- nama
  vale <- which( skel$G == 1 )
  dag <- NULL
  score <- NULL
  if ( length(vale) > 0 ) {
    x <- as.data.frame(x)
    mhvale <- as.data.frame( which(skel$G == 0, arr.ind = TRUE) )
    mhvale[, 1] <- nama[ mhvale[, 1] ]
    mhvale[, 2] <- nama[ mhvale[, 2] ]
    colnames(mhvale) <- c("from", "to")
    if ( !is.null(blacklist) ) {
      colnames(blacklist) <- c("from", "to")
      mhvale <- rbind(mhvale, blacklist)
    }
    if ( !is.null(whitelist) )  colnames(whitelist) <- c("from", "to")
    if ( algo == "hc" ) {
      dag <- bnlearn::hc(x, blacklist = mhvale, whitelist = whitelist, score = score, restart = restart)
    } else {
      dag <- bnlearn::tabu(x, blacklist = mhvale, whitelist = whitelist, score = score, tabu = tabu, restart = restart)
    }
    score <- bnlearn::score(x = dag, data = x, type = score )
  }
  runtime <- proc.time() - runtime
  list(ini = skel, dag = dag, score = score, runtime = runtime)
}












