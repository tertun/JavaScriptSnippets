\name{gen.ker}
\alias{gen.ker}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generate Kernel Matrix
}
\description{
Generate kernel matrices three commonly used kernel functions in practice, including Gaussian, linear and exponential kernels. Output is n by n kernel matrix without centralization.
}
\usage{
gen.ker(covx, rho, kernel.index = c("e","l","g"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{covx}{
input matrix, of dimension n by p; each row is an p-dim observation vector.
}
  \item{rho}{
tuning parameter in kernel functions.
}
\item{kernel.index}{
kernel function index, where "e" indicates exponential kernel, "l" indicates linear kernel, and "g" indicates gaussian kernel. See explicit function expressions in details section.
}
}
\details{
The three candidate kernel functions are defined as follows, where \eqn{p} is the dimension of the covariants and \eqn{rho} is tuning parameter

Linear kernel: \deqn{K_L(x1,x2)=x1^T x2/rho}
Gaussion kernel: \deqn{K_G(x1,x2)=exp(-||x1-x2||^2/rho)}
Exponential kernel: \deqn{K_E(x1,x2)=exp{-(||x1||^2+3||x1-x2||^2+||x2||^2)/rho}}
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
\examples{
n <- 100
p <- 500
covx <- matrix(rnorm(n*p),n,p)
beta <- c(rep(0,480),rnorm(20,0,0.5))
y <- rnorm(n)+covx\%*\%beta
ker.e <- gen.ker(covx,p,"e")
ker.l <- gen.ker(covx,p,"l")
ker.g <- gen.ker(covx,p,"g")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
