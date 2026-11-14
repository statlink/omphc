omp.network.boot <- function(x, R = NULL, method = "pvalue", tol = 0.05, B = 200) {
  G <- omphc::omp.network(x = x, R = R, method = method, tol = tol)$G
  dm <- dim(x)
  n <- dm[1]   ;     p <- dm[2]

  runtime <- proc.time()
  Gboot <- matrix(0, p, p)
  for (i in 1:B) {
    id <- Rfast2::Sample.int(n, n, replace = TRUE)
    Gboot <- Gboot + omphc::omp.network(x = x[id, ], method = method, tol = tol)$G
  }  ## end for (i in 1:B)
  runtime <- proc.time() - runtime

  list(G = G, Gboot = Gboot/B, runtime = runtime)
}
