\name{KerRegSel}
\alias{KerRegSel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Selection of Kernel Regularization Parameter
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
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
\item{ker}{n by n kernel matrix, noncentralized or centralized. It can be self defined or the ouput of built-in function gen.ker().
}
\item{lambda.seq}{A (positive) sequence of regularization parameter lambda.}

\item{L}{Number of folds that samples divided into in the proposed stability selection procedure. }
}

\details{
%%  ~~ If necessary, more details than the description above ~~
}


\value{
\item{lambda.stab.sel}{the tuning parameter selected in the stability procedure.}
\item{lambda.sel.list}{the complete list of selected lambdas when the \emph{l}th fold was left out, where \eqn{l=1,2,...,L}.
}
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
ker <- gen.ker(covx,p,"l")

KerRegSel(y,covx,ker)
KerRegSel(y,covx,ker,lambda.seq=10^seq(-3,5,by=0.5),L=10)
}
