\name{ABOmle}
\alias{ABOmle}
\title{
Compute the MLE of gene frequencies of A,B,O
}
\description{
Given observed counts nO,nA,nB,nAB,
compute the maximum likelihood estimates of p, q, r
(corresponding to the gene frequencies of A,B,O),
using the Fisher scoring or the Newton-Raphson methods.
}
\usage{
ABOmle(nO, nA, nB, nAB, p = 0.3, q = 0.3, acu = 1e-3, max_iter = 10, method = "N", choose = FALSE)
}
\arguments{
  \item{nO,nA,nB,nAB}{
Observed counts of phenotype O,A,B,AB
}
  \item{p, q}{
The initial value of the iteration
}
  \item{acu}{
Iteration accuracy
}
  \item{max_iter}{
The maximal number of iterations
}
  \item{method}{
The method to be used;
for Newton-Raphson method, input value "N";
for Fisher-Scoring, input value "F";
}
  \item{special}{
whether or not choose special p, q to have better convergence speed and stability
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  \item{r}{The maximum likelihood estimates of r(the gene frequencies of O)}
  \item{p}{The maximum likelihood estimates of p(the gene frequencies of A)}
  \item{q}{The maximum likelihood estimates of q(the gene frequencies of B)}
  \item{iter}{The number of iterations}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Tongyue Jia, Qilin Li, Yejing Wang
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
  ABOmle(153,234,456,387, method = "N")
}



% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
