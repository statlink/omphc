update.skel <- function(skel, method = "pvalue", tol = 0.01) {

  p <- length(skel$Rho)
  G <- matrix(0, p, p)
  colnames(G) <- rownames(G) <- colnames(skel$G)

  if ( method == "pvalue" ) {
    tol <- log(tol)
    for ( i in 1:p ) {
      a <- skel$Rho[[ i ]][-1]
      id <- a < tol
      b <- match(FALSE, id) - 1
      if ( !is.na(b) ) {
        sel <- as.numeric( names(id) )[1:b]
        skel$Rho[[ i ]][-1] <- a[1:b]
      } else  sel <- as.numeric( names(id) )
      G[i, sel] <- 1
    }

  } else {
    for ( i in 1:p ) {
      a <- skel$Rho[[ i ]]
      id <- abs( diff(a) ) > tol
      b <- match(FALSE, id)
      if ( !is.na(b) ) {
        sel <- as.numeric( names(id) )[1:c(b-1)]
        skel$Rho[[ i ]] <- a[1:b]
      } else  sel <- as.numeric( names(id) )
      G[i, sel] <- 1
    }
  }
  a <- which( G == 1  &  t(G) == 1 )
  G[ -a ] <- 0
  skel$G <- G
  skel
}
