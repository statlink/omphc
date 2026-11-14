\name{omphc}
\alias{omphc}
\title{
The OMP-based Bayesian network learning algorithm
}

\description{
The OMP-based Bayesian network learning algorithm.
}

\usage{
omphc(x, skel = NULL, R = NULL, method = "pvalue", tol = 0.05, algo = "hc",
tabu = 10, restart = 10, score = "bic-g", blacklist = NULL, whitelist = NULL)
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
The scoring algorithm, either "hc" or "tabu"
}
\item{tabu}{
An integer, the length of the tabu list used in the tabu function.
}
\item{restart}{
An integer, the number of random restarts.
}
\item{score}{
A character string, the label of the network score to be used in the algorithm. If none is specified,
the default score is the Bayesian Information Criterion for continuous data sets.
The available score for continuous variables are: "bic-g" (default), "loglik-g", "aic-g", "bic-g" or "bge".
}
\item{blacklist}{
A data frame with two columns (optionally labeled "from" and "to"), containing a set of arcs not to
be included in the graph.
}
\item{whitelist}{
A data frame with two columns (optionally labeled "from" and "to"), containing a set of arcs to be
included in the graph.
}
}

\details{
The OMPHC or OMPTABU algorithm is implemented. The OMP algorithm (Tsagris et al., 2002,
Mallat S. G. and Zhang Z. (1993)), is implemented during the skeleton identification phase.
Next, the Hill Climbing greedy search or the Tabu search is employed to score the network.
}

\value{
A list including:
\item{ini}{
A list including the output of the \code{ \link{omp.network} } function.
}
\item{dag}{
A "bn" class output. A list including the outcome of the Hill-Climbing or the Tabu search phase.
See the package "bnlearn" for more details.
}
\item{scoring}{
The score value.
}
\item{runtime}{
The duration of the algorithm.
}
}

\references{
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
\code{ \link{omp.network}
}
}

\examples{
m <- rbn2(1000, p = 50, nei = 5)
G <- m$G
x <- m$x
a <- omphc(x, tol = 0.01)
}
