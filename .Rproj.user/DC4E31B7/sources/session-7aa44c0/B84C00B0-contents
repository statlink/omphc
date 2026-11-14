\name{fast.hc}
\alias{fast.hc}
\title{
Fast hill climbing for Gaussian Bayesian networks with the BIC score
}

\description{
Greedy hill-climbing structure learning of a Gaussian Bayesian network using
the BIC score, computed entirely from sufficient statistics (the correlation
matrix). The search uses the same moves (arc addition, deletion and reversal),
the same move order, tie-breaking and tolerance as \code{bnlearn::hc()}, and it returns the same network as
\code{bnlearn::hc(x, score = "bic-g")} (i.e. with no random restarts).
}

\usage{
fast.hc(x, R, blacklist = NULL, whitelist = NULL)
}

\arguments{
\item{x}{
A numeric matrix with continuous variables and no missing values.
}
\item{R}{
The correlation matrix of \code{x}, e.g. \code{cor(x)}, with its columns in the
same order as those of \code{x}.
}
\item{blacklist}{
A two-column matrix of arcs (from, to) that must not be in the
network, or \code{NULL}.
}
\item{whitelist}{
A two-column matrix of arcs (from, to) that must be in the
network, or \code{NULL}. The search starts from these arcs, they are never
removed, and their reverses are blacklisted.
}
}

\details{
The data are used only once, to compute the \eqn{p \times p} correlation
matrix \eqn{\mathbf{R}}. For a node \eqn{j} with parent set \eqn{P} of size
\eqn{k}, let \eqn{\mathbf{R}_{PP} = \mathbf{L}\mathbf{L}^\top} and
\eqn{\mathbf{B} = \mathbf{L}^{-1}\mathbf{R}_{P\cdot}}. The residual sum of squares of
every model obtained by adding a candidate parent \eqn{i} follows from the
partial covariances and variances,
\eqn{RSS_j(P \cup \{i\}) = RSS_j(P) - c_i^2/d_i}, with
\eqn{c_i = R_{ij} - \mathbf{b}_i^\top\mathbf{b}_j} and
\eqn{d_i = 1 - \mathbf{b}_i^\top\mathbf{b}_i}. Deleting a parent \eqn{m} gives
\eqn{RSS_j(P \setminus \{m\}) = RSS_j(P) + \beta_m^2 / (\mathbf{R}_{PP}^{-1})_{mm}}.
After each move, the cached score differences of the affected node(s) are
recomputed in \eqn{O(k^3 + pk^2)} operations, independently of the sample
size. Acyclicity is checked lazily, by a depth-first search, only for
candidate moves that improve on the current best.
}
\value{
A list with components:
\item{arcs}{A two-column character matrix with the arcs of the learned network.}
\item{amat}{The adjacency matrix of the network.}
\item{score}{The BIC score of the network (on the log-likelihood scale, i.e.
larger is better, as in \code{bnlearn}).}
\item{node.score}{The BIC score of each node.}
\item{iterations}{The number of moves performed.}
\item{moves}{The number of additions, deletions and reversals performed.}
\item{bn}{The network as a \code{bn} object of the \pkg{bnlearn} package.}
\item{runtime}{The run time of the function.}
}
\references{
Scutari M. (2010). Learning Bayesian networks with the bnlearn R package.
Journal of Statistical Software, 35(3): 1--22.

Tsamardinos I., Brown L.E. and Aliferis C.F. (2006). The max-min hill-climbing
Bayesian network structure learning algorithm. Machine Learning, 65(1): 31--78.
}

\author{
Michail Tsagris.
}
\seealso{
\code{\link{omphc}, \link{omp.network}}
}

\examples{
x <- matrix(rnorm(100 * 100), ncol = 100)
mod <- fast.hc(x, R = cor(x))
mod$arcs
mod$score
}

