\name{baymode}
\alias{baymode}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Bayesian multivariate mode estimation
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
baymode(X,niter=1000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
%%     ~~Describe \code{x} here~~
Data set, a n x d matrix.
}
\item{niter}{
Number of iterations. Default: 1000.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{mu}{niter simulations of the multivariate mode.}
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
\keyword{ ~mode }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~bayesian }% __ONLY ONE__ keyword per line
