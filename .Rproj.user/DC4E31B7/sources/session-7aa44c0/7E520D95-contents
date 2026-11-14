omphc.fast <- function(x, skel = NULL, R = NULL, method = "pvalue", tol = 0.05, blacklist = NULL, whitelist = NULL) {
  runtime <- proc.time()
  if ( is.null(skel) )  skel <- omphc::omp.network(x, R = R, method = method, tol = tol)
  nama <- colnames(skel$G)
  if ( is.null(nama) )  nama <- paste("X", 1:dim(x)[2], sep = "")
  colnames(x) <- nama
  vale <- which( skel$G == 1 )
  R <- skel$R
  diag(R) <- 1
  dag <- NULL
  score <- NULL
  if ( length(vale) > 0 ) {
    mhvale <- as.data.frame( which(skel$G == 0, arr.ind = TRUE) )
    b <- which( mhvale[, 1] - mhvale[, 2] == 0)
    mhvale <- mhvale[-b, ]
    mhvale[, 1] <- nama[ mhvale[, 1] ]
    mhvale[, 2] <- nama[ mhvale[, 2] ]
    colnames(mhvale) <- c("from", "to")
    if ( !is.null(blacklist) ) {
      colnames(blacklist) <- c("from", "to")
      mhvale <- rbind(mhvale, blacklist)
    }
    if ( !is.null(whitelist) )  colnames(whitelist) <- c("from", "to")
    dag <- omphc::fast.hc(x, R = R, blacklist = mhvale, whitelist = whitelist)
  }
  runtime <- proc.time() - runtime
  list(ini = skel, dag = dag, runtime = runtime)
}
