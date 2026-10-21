\name{update.skel}
\alias{update.skel}
\title{
Update the skeleton of the Bayesian network
}

\description{
The OMP-based Bayesian network learning algorithm.
}

\usage{
update.skel(skel, method = "pvalue", tol = 0.01)
}

\arguments{
\item{skel}{
The output of the function \code{\link{op.network}}.
}
\item{method}{
The method used to construct the skeleton, p-value based ("pvalue"), BIC based ("BIC"),
or adjusted R-square ("ar2") based.
}
\item{tol}{
The tolerance value to update the skeleton. The significance level, a value for the BIC difference,
or a percentage for the adjusted R-square. It has to be smaller than the tolerance value initially used.
}
}

\details{
Suppose you have run the skeleton of the OMPHC with a tolerance value equal to 0.05 using the "pvalue"
method. Now, you wish to update the skeleton for tol=0.01. Instead of running the \code{ \link{omp.network} }
again with tol=0.01, you can update the G matrix with this function. The function will change the G
matrix, and the Rho list.
}

\value{
The same output as \code{ \link{omp.network} }.
}

\references{
Tsagris M., Papadovasilakis Z., Lakiotaki K. and Tsamardinos I. (2022). The \eqn{\gamma}-OMP
Algorithm for Feature Selection With Application to Gene Expression Data.
IEEE/ACM Transactions on Computational Biology and Bioinformatics, 19(2): 1214--1224.

Tsagris M. (2022). The FEDHC Bayesian Network Learning Algorithm. Mathematics, 10(15): 2604.

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
a1 <- omp.network(x, tol = 0.1)
a2 <- omp.network(x, tol = 0.01)
b <- update.skel(a1, method = "pvalue", tol = 0.01)
}
