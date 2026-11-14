omphc.boot <- function(x, skel = NULL, R = NULL, method = "pvalue", tol = 0.05, algo = "hc",
                       tabu = 10, restart = 10, score = "bic-g", blacklist = NULL, whitelist = NULL,
                       B = 200, ncores = 1) {

  mod <- omphc::omphc(x, skel = skel, R = R, method = method, tol = tol, algo = algo,
                      tabu = tabu, restart = restart, score = score,
                      blacklist = blacklist, whitelist = whitelist)
  dm <- dim(x)
  n <- dm[1]   ;     p <- dm[2]

  ## one bootstrap replicate, returns the vectorised adjacency matrix
  boot1 <- function(i) {
    id <- Rfast2::Sample.int(n, n, replace = TRUE)
    gb <- omphc::omphc(x[id, ], method = method, tol = tol, algo = algo, tabu = tabu, restart = restart,
                       score = score, blacklist = blacklist, whitelist = whitelist)
    as.vector( omphc::bnmat(gb$dag) )
  }

  runtime <- proc.time()
  if ( ncores == 1 ) {
    Gboot <- numeric(p * p)
    for ( i in 1:B )  Gboot <- Gboot + boot1(i)
  } else {
    cl <- parallel::makePSOCKcluster(ncores)
    on.exit( parallel::stopCluster(cl), add = TRUE )
    Gboot <- parallel::parSapply(cl, 1:B, boot1)   ## (p^2) x B matrix
    Gboot <- Rfast::rowsums(Gboot)
  }
  runtime <- proc.time() - runtime
  Gboot <- matrix(Gboot, nrow = p, ncol = p)

  list(mod = mod, Gboot = Gboot / B, runtime = runtime)
}
