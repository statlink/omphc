fast.hc <- function(x, R, blacklist = NULL, whitelist = NULL) {

  runtime <- proc.time()

  n <- dim(x)[1]   ;   p <- dim(x)[2]
  nodes <- colnames(x)
  if (is.null(nodes)) nodes <- paste0("X", seq_len(p))

  xc <- x - rep( Rfast::colmeans(x), each = n )
  css <- Rfast::colsums(xc^2)
  rm(xc)

  ## arc sets -> adjacency matrices.
  arcs2amat <- function(a, what) {
    m <- matrix(0L, p, p)
    if (is.null(a) || NROW(a) == 0) return(m)
    a <- as.matrix(a)
    if (ncol(a) != 2) stop(what, " must be a two-column matrix.")
    idx <- matrix(match(as.character(a), nodes), ncol = 2)
    if (anyNA(idx)) stop(what, " contains unknown node names.")
    if (any(idx[, 1] == idx[, 2])) stop(what, " contains loops.")
    m[idx] <- 1L
    m
  }
  bl <- arcs2amat(blacklist, "blacklist")
  wl <- arcs2amat(whitelist, "whitelist")
  if (any(bl & wl)) stop("some arcs are both whitelisted and blacklisted.")
  if (any(wl & t(wl))) stop("the whitelist contains arcs in both directions.")
  ## as in bnlearn, the reverse of a whitelisted arc is blacklisted.
  bl <- bl | t(wl)
  storage.mode(bl) <- "integer"

  ## last argument FALSE = unbiased variance, i.e. bnlearn's "bic-g".
  maxp <- Inf   ;  max.iter = Inf
  res <- .fasthc_cpp(R, css, n, bl, wl, Inf, Inf, FALSE)

  ## node and network BIC scores (bnlearn's "bic-g" scale).
  k <- res$nparents
  rss <- res$rss * css
  df <- n - k - 1
  loglik <- -0.5 * n * log(2 * pi * rss / df) - 0.5 * df
  node.score <- loglik - 0.5 * log(n) * (k + 2)
  names(node.score) <- nodes

  amat <- res$amat
  dimnames(amat) <- list(nodes, nodes)
  w <- which(amat == 1L, arr.ind = TRUE)
  w <- w[order(w[, 1], w[, 2]), , drop = FALSE]
  arcs <- cbind(from = nodes[w[, 1]], to = nodes[w[, 2]])

  out <- list(arcs = arcs, amat = amat, score = sum(node.score),
              node.score = node.score, iterations = res$iterations,
              moves = structure(res$moves, names = c("add", "drop", "reverse")))

  g <- bnlearn::empty.graph(nodes)
  bnlearn::arcs(g) <- arcs
  out$bn <- g

  out$runtime <- proc.time() - runtime
  out

}
