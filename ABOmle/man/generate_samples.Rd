\name{generate_samples}
\alias{generate_samples}
\title{
Generate observed counts of four blood groups
}
\description{
Given p, q(the gene frequencies of A, B),
generate observed counts of four blood groups
(O, A, B, AB)
}
\usage{
generate_samples(n = 1e+04, p = 0.3, q = 0.3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
The size of sample
}
  \item{p, q}{
the gene frequencies of A, B
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
return a table with observed counts of four blood groups
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
generate(5000, 0.6 0.2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
