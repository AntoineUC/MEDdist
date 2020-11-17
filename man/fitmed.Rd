\name{fitmed}
\alias{fitmed}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Maximum Likelihood estimation of the parameters of a MED distribution
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function estimates the parameters mu, Sigma, rho and nu of the MED distribution.
}
\usage{
fitmed(x,mu0,Sigma0,rho0,nu0,maxit)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{X} here~~
Data. A n x p matrix, p being the dimension, n the size of the sample.
}
  \item{mu0}{
%%     ~~Describe \code{k} here~~
p-dimensional vector. Initial value of mu. Default: Empirical mean of x.
}
\item{Sigma0}{
%%     ~~Describe \code{k} here~~
p x p matrix. Initial value of Sigma. Default: Empirical variance of x.
}
\item{rho0}{
%%     ~~Describe \code{k} here~~
Real value between 0 and 1. Intial value of the asymmetry parameter rho. Default: 0.1.
}
\item{nu0}{
%%     ~~Describe \code{k} here~~
p-dimensional vector. Initial value of the asymmetry direction nu. Default: e1.
}
\item{maxit}{
%%     ~~Describe \code{k} here~~
Maximum number of iterations. Default: 100.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
A list of the estimated parameters.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
Antoine Usseglio-Carleve
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
\keyword{ ~MED }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~Maximum Likelihood }% __ONLY ONE__ keyword per line
