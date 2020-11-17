\name{dmed}
\alias{dmed}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Probability density function of the MED distribution
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function computes the probability density function of the MED distribution.
}
\usage{
dmed(x,mu,Sigma,rho,nu)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{X} here~~
p x n matrix, p being the dimension, n the number of points where the p.d.f. is evaluated.
}
  \item{mu}{
%%     ~~Describe \code{k} here~~
p-dimensional vector of expectile. Default: rep(0,p).
}
\item{Sigma}{
%%     ~~Describe \code{k} here~~
p x p matrix. Default: diag(p).
}
\item{rho}{
%%     ~~Describe \code{k} here~~
Asymmetry parameter. Default: 0.
}
\item{nu}{
%%     ~~Describe \code{k} here~~
Asymmetry direction. Default: rep(1,p)/sqrt(p).
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
A vector of the p.d.f. evaluated at n points.
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
\keyword{ ~density }% __ONLY ONE__ keyword per line
