\name{rmed}
\alias{rmed}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Simulation of the MED distribution
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function returns a sample of the MED distribution.
}
\usage{
rmed(n,mu,Sigma,rho,nu)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
%%     ~~Describe \code{X} here~~
Size of the sample.
}
  \item{mu}{
%%     ~~Describe \code{k} here~~
p-dimensional vector of expectile. Default: 0.
}
\item{Sigma}{
%%     ~~Describe \code{k} here~~
p x p matrix. Default: diag(length(mu)).
}
\item{rho}{
%%     ~~Describe \code{k} here~~
Asymmetry parameter. Default: 0.
}
\item{nu}{
%%     ~~Describe \code{k} here~~
Asymmetry direction. Default: e1.
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
A sample of size n of the MED distribution.
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
\keyword{ ~simulation }% __ONLY ONE__ keyword per line
