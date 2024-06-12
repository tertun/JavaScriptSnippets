\name{KerUTest}
\alias{KerUTest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Testing Nonparametric Functions in RKHS %%  ~~function to do ... ~~
}
\description{
Testing nonparametric functions in Reproducing Kernel Hilbert Space (RKHS) of high dimensional or functional covariates. A U-statistics was developed to test the relationship between the quntitative response variable and a p-dimensional covariates.
}
\usage{
KerUTest(y, covx, ker, alpha=0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
p-dim quantitative reponse vector.
}
\item{covx}{
input matrix, of dimension n by p; each row is an p-dim observation vector.
}
%\item{kernel.index}{
%kernel function index, where "e" indicates exponential kernel, "l" indicates linear kernel, and "g" indicates gaussian kernel. See explicit function expressions in details section.
%}
\item{ker}{n by n kernel matrix, noncentralized or centralized. It can be self defined or the ouput of built-in function gen.ker().
}
\item{alpha}{
this argument is used to estimate the testing power of alpha-level test. Default setting is alpha=0.05.
}
}
\details{Scaled chi-square approximation works for any type of dimension (low or high). Normal approximation works better in high dimensional setting than low dimensional setting. Kernel functions are defined as follows:
Linear kernel: \deqn{K_L(x1,x2)=x1^T x2/p}
Gaussion kernel: \deqn{K_G(x1,x2)=exp(-||x1-x2||^2/p)}
Exponential kernel: \deqn{K_E(x1,x2)=exp{-(||x1||^2+3||x1-x2||^2+||x2||^2)/p}
}}


\value{
\item{pv.norm}{p-values calcluated based on normal distribution.}
\item{pv.chisq}{p-values calcluated based on scaled chi-square distribution.}
\item{power.chisq}{estimated power baed on scaled chi-square distribution.}
\item{shift}{the signal-to-noise-ratio (SNR). The larger shift indicates stronger signal of association.}
}
\references{
Tao He, Shaoyu Li, Ping-Shou Zhong and Yuehua Cui. "An optimal kernel‐based U‐statistic method for quantitative gene‐set association analysis." \emph{Genetic Epidemiology} 43.2 (2019): 137-149.

Tao He, Ping-Shou Zhong, Yuehua Cui and Vidyadhar Mandrekar. "Unified Tests for Nonparametric Functions in RKHS with Kernel Selection and Regularization". \emph{Submitted.}
}
\author{
Tao He, Ping-Shou Zhong, Yuehua Cui and Vidyadhar Mandrekar
Maintainer: Tao He \email{hetao@sfsu.edu}
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{n <- 100
p <- 500
covx <- matrix(rnorm(n*p),n,p)
beta <- c(rep(0,480),rnorm(20,0.5,1))
y <- rnorm(n)+covx\%*\%beta

#testing using exponential kernel
ker.e <- gen.ker(covx,p,"e")
test.e <- kerUTest(y,covx,ker.e)
test.e

#testing using linear kernel
ker.l <- gen.ker(covx,p,"l")
test.l <- kerUTest(y,covx,ker.l)
test.l

#testing using Gaussian kernel
ker.g <- gen.ker(covx,p,"g")
test.g <- kerUTest(y,covx,ker.g)
test.g}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
