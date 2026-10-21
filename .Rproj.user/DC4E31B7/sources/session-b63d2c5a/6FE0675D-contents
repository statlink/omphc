\name{omp.network}
\alias{omp.network}
\title{
skeleton of a Bayesian network produced by the OMP algorithm
}

\description{
skeleton of a Bayesian network produced by the OMP algorithm.
}

\usage{
omp.network(x, R = NULL, method = "pvalue", tol = 0.05)
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
}

\details{
As in FEDHC the first phase consists of a variable selection procedure, the OMP algorithm.
Then, hill-climbing or TABU search is applied to direct the edges.
}

\value{
A list including:
\item{ini.stat}{
The test statistics of the univariate associations.
}
\item{ini.pvalue}{
The initial p-values univariate associations.
}
\item{pvalue}{
A matrix with the logarithm of the p-values of the updated associations.
This final p-value is the maximum p-value among the two p-values in the end.
}
\item{runtime}{
The duration of the algorithm.
}
\item{ntests}{
The number of tests conducted during each k.
}
\item{G}{
The adjancency matrix. A value of 1 in G[i, j] appears in G[j, i] also, indicating that i and j
have an edge between them.
}
}
\references{
Tsagris M. (2022). The FEDHC Bayesian Network Learning Algorithm.
Mathematics, 10(25): 2604.

Szekely G.J., Rizzo M.L. and Bakirov N.K. (2007). Measuring and Testing Independence
by Correlation of Distances. Annals of Statistics, 35(6): 2769--2794.

Szekely G.J. and Rizzo M. L. (2014). Partial distance correlation with methods for dissimilarities.
Annals of Statistics, 42(6): 2382--2412.

Huo X. and Szekely G.J. (2016). Fast computing for distance covariance.
Technometrics, 58(4): 435--447.

Tsamardinos I., Brown E.L. and Aliferis F.C. (2006). The max-min hill-climbing Bayesian network
structure learning algorithm. Machine Learning, 65(1): 31--78.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}
\seealso{
\code{\link{omphc}
}
}

\examples{
m <- rbn2(1000, p = 50, nei = 5)
G <- m$G
x <- m$x
a <- omp.network(x, tol = 0.1)
}

