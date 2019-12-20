\name{generate_samples_simulation}
\alias{generate_samples_simulation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generate observed counts of four blood groups
in the way how genes interact with each other
}
\description{
Given p, q(the gene frequencies of A, B),
generate observed counts of four blood groups
(O, A, B, AB)
in the way how genes interact with each other
}
\usage{
generate_samples_simulation (size = 100,total = 2, W=c("A","B","O"),wt=c(1,1,0),wp=c(0.8,0.1,0.1))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{size}{
The size of sample
}
  \item{total}{
the total number of attributes in each individual
}
  \item{w}{
all the possibilities of the attributes
}
  \item{wt}{
the corresponding weight
}
  \item{wp}{
the corresponding frequency
}
}
\details{
This is function simuating the process of generating blood groups, and produces result indifferent from 'generate.samples', but it is more close to the real case. 'size' represents the size of the sample. 'total' is the total number of attributes in each individual. 'W' contains all the possibilities of the attributes. 'wt' and 'wp' are the corresponding weight and frequency. When total, W and wt are the default arguments, this function is identical to the simulation of blood groups. This function may apply to other occasions as well.
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
generate.samples_simulation()
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
