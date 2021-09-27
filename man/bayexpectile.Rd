\name{bayexpectile}
\alias{bayexpectile}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Bayesian multivariate expectile estimation
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
bayexpectile(X,dir=(1,0,...,0),rho=0.5,niter=100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
%%     ~~Describe \code{x} here~~
Data set, a n x d matrix.
}
\item{dir}{
Direction of the multivariate expectile level. Must be on the unit circle. Default: (1,0,...,0).
}
\item{rho}{
Norm of the multivariate expectile. Default: 0.5.
}
\item{niter}{
Number of iterations. Default: 100.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{mu}{niter simulations of the multivariate expectile.}
  \item{muhat}{Mean of the niter simulations.}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~expectile }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~bayesian }% __ONLY ONE__ keyword per line
