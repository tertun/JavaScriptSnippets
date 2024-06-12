\name{KerSel}
\alias{KerSel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Kernel Selection
}
\description{
Among Gaussian, linear and exponential kernels, choose the one that has best signal-to-noise-ratio (SNR) so that the proposed nonparametric testing is most powerful.
}
\usage{
KerRegSel(y, covx, ker, lambda.seq=10^seq(-5,5,by=1), L=5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
p-dim quantitative reponse vector.
}
\item{covx}{
input matrix, of dimension n by p; each row is an p-dim observation vector.
}
}
\details{The three candidate kernel functions are defined as follows, where \eqn{p} is the dimension of the covariants. The kernel with the largest SNR is the selected kernel function.

Linear kernel: \deqn{K_L(x1,x2)=x1^T x2/p}
Gaussion kernel: \deqn{K_G(x1,x2)=exp(-||x1-x2||^2/p)}
Exponential kernel: \deqn{K_E(x1,x2)=exp{-(||x1||^2+3||x1-x2||^2+||x2||^2)/p}
}}
\value{
\item{shift.e}{the signal-to-noise-ratio (SNR) if testing using exponential kernel.}
\item{shift.l}{the signal-to-noise-ratio (SNR) if testing using linear kernel.}
\item{shift.g}{the signal-to-noise-ratio (SNR) if testing using Gaussian kernel.}
}
\references{
Tao He, Shaoyu Li, Ping-Shou Zhong and Yuehua Cui. "An optimal kernel‐based U‐statistic method for quantitative gene‐set association analysis." \emph{Genetic Epidemiology} 43.2 (2019): 137-149.

Tao He, Ping-Shou Zhong, Yuehua Cui and Vidyadhar Mandrekar. "Unified Tests for Nonparametric Functions in RKHS with Kernel Selection and Regularization". \emph{Submitted.}
}
\author{
Tao He, Ping-Shou Zhong, Yuehua Cui and Vidyadhar Mandrekar
Maintainer: Tao He \email{hetao@sfsu.edu}
}

\examples{
n <- 100
p <- 500
covx <- matrix(rnorm(n*p),n,p)
beta <- c(rep(0,480),rnorm(20,0,0.5))
y <- rnorm(n)+covx\%*\%beta
KerSel(y,covx)
}

