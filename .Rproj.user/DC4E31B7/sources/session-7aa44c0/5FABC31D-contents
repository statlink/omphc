\name{omp.network.boot}
\alias{omp.network.boot}
\title{
Bootstrap versions of the skeleton of the OMPHC
}

\description{
Bootstrap versions of the skeleton of the OMPHC.
}

\usage{
omp.network.boot(x, R = NULL, method = "pvalue", tol = 0.05, B = 200)
}

\arguments{
\item{x}{
A numerical matrix with the variables. If you have a data.frame (i.e. categorical data) turn them
into a matrix. Note, that for the categorical case data,
the numbers must start from 0. No missing data are allowed.
}
\item{R}{
If the correlation matrix is available, pass it here.
}
\item{method}{
The method to proceed with the OMP, p-value based ("pvalue"), BIC based ("BIC"), or adjusted R-square
("ar2") based.
}
\item{tol}{
The tolerance value to terminate the OMP. The significance level, a value for the BIC difference,
or a percentage for the adjusted R-square.
}
\item{B}{
The number of bootstrap resamples to draw. The algorithm is performed in each bootstrap sample.
In the end, the adjacency matrix on the observed data is returned, along with another adjacency
matrix produced by the bootstrap. The latter one contains values from 0 to 1 indicating the
proportion of times an edge between two nodes was present.
}
}

%details{
%
%}

\value{
A list including:
\item{G}{
The observed adjancency matrix. A value of 1 in G[i, j] appears in G[j, i] also, indicating that
i and j have an edge between them.
}
\item{Gboot}{
The bootstrapped adjancency matrix. A value of 1 in G[i, j] appears in G[j, i] also, indicating
that i and j have an edge between them.
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
\code{ \link{omp.network}, \link{omphc}
}
}

\examples{
x <- rbn2(500, p = 20, nei = 3)$x
a <- omp.network.boot(x)
}
