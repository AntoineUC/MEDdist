\name{qmed}
\alias{qmed}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Quantiles of the MED distribution
  %%  ~~function to do ... ~~
}
\description{
  %%  ~~ A concise (1-5 lines) description of what the function does. ~~
    This function computes the quantiles of the MED distribution.
}
\usage{
  pmed(probs,mu,Sigma,rho,nu)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
    %%     ~~Describe \code{X} here~~
      vector of probabilities.
  }
  \item{mu}{
    %%     ~~Describe \code{k} here~~
      Value of expectile. Default: 0.
  }
  \item{Sigma}{
    %%     ~~Describe \code{k} here~~
      Scale parameter Sigma. Default: 1.
  }
  \item{rho}{
    %%     ~~Describe \code{k} here~~
      Asymmetry parameter. Default: 0.
  }
  \item{nu}{
    %%     ~~Describe \code{k} here~~
      Asymmetry direction. Default: 1.
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
  A vector of quantiles.
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
\keyword{ ~quantile }% __ONLY ONE__ keyword per line