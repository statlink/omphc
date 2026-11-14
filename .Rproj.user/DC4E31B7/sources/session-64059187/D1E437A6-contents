\name{omphc.boot}
\alias{omphc.boot}
\title{
Bootstrapping the OMPHC Bayesian network learning algorithm
}

\description{
Bootstrapping the OMPHC Bayesian network learning algorithm.
}

\usage{
omphc.boot(x, skel = NULL, R = NULL, method = "pvalue", tol = 0.05, algo = "hc",
tabu = 10, restart = 10, score = "bic-g", blacklist = NULL, whitelist = NULL,
B = 200, ncores = 1)
}

\arguments{
\item{x}{
A numerical matrix with the variables. No missing data are allowed.
}
\item{skel}{
If you have the output of the skeleton phase, the output from the function
\code{\link{omp.network}} plug it here.
This can save time.
}
\item{R}{
If the correlation matrix is available, pass it here.
}
\item{method}{
If the skeleton is not available, the method to proceed with the OMP,
p-value based ("pvalue"), BIC based ("BIC"), or adjusted R-square ("ar2") based.
}
\item{tol}{
The tolerance value to terminate the OMP. The significance level, a value for the BIC difference,
or a percentage for the adjusted R-square.
}
\item{algo}{
The scoring algorithm, either "hc" or "tabu".
}
\item{tabu}{
An integer, the length of the tabu list used in the tabu function.
}
\item{restart}{
An integer, the number of random restarts.
}
\item{score}{
A character string, the label of the network score to be used in the algorithm.
If none is specified, the default score is the Bayesian Information Criterion for
continuous data sets. The available score for continuous variables are:
"bic-g" (default), "loglik-g", "aic-g", "bic-g" or "bge".
}
\item{blacklist}{
A data frame with two columns (optionally labeled "from" and "to"), containing
a set of arcs not to be included in the graph.
}
\item{whitelist}{
A data frame with two columns (optionally labeled "from" and "to"), containing
a set of arcs to be included in the graph.
}
\item{B}{
The number of bootstrap resamples to draw. The algorithm is performed in each
bootstrap sample. In the end, the adjacency matrix on the observed data is returned,
along with another adjacency matrix produced by the bootstrap. The latter one
contains values from 0 to 1 indicating the proportion of times an edge between
two nodes was present.
}
\item{ncores}{
The number of cores to use, in case of parallel computing.
}
}

\details{
The OMPHC algorithm is implemented for every bootstrap sample.
}

\value{
A list including:
\item{mod}{
A list including the output of the \code{ \link{omphc} } function.
}
\item{Gboot}{
The bootstrapped adjancency matrix of the Bayesian network.
}
\item{runtime}{
The duration of the algorithm.
}
}

\references{
Tsagris M. and Adam C. (2026). The OMPHC Bayesian network learning algorithm.
https://www.researchgate.net/profile/Michail-Tsagris/publication/416700162_The_OMPHC_Bayesian_network_learning_algorithm/links/6adb9a3b21cbf015d2a094a6/The-OMPHC-Bayesian-network-learning-algorithm.pdf

Tsagris M., Papadovasilakis Z., Lakiotaki K. and Tsamardinos I. (2022). The \eqn{\gamma}-OMP
Algorithm for Feature Selection With Application to Gene Expression Data.
IEEE/ACM Transactions on Computational Biology and Bioinformatics, 19(2): 1214--1224.

Tsamardinos I., Brown E.L. and Aliferis F.C. (2006). The max-min hill-climbing Bayesian network
structure learning algorithm. Machine Learning, 65(1): 31--78.

Mallat S. G. and Zhang Z. (1993). Matching pursuits with time-frequency dictionaries.
IEEE Transactions on Signal Processing, 41(12): 3397--3415.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}
\seealso{
\code{ \link{omphc} }
}

\examples{
x <- matrix( rnorm(200 * 20, 1, 10), nrow = 200 )
a <- omphc.boot(x, B = 50)
}
