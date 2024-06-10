#remove all the objects
rm(list=ls(all=TRUE))

#checking packages
list.of.packages <- c("psych","caret")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

#loading packages
#require("mvtnorm")
#require("ttutils")
#require("methods")
library("psych")
library("caret")
#require("Runuran")
#require("EQL")

#define kernel functions
##########function to generate kernel matrix###############
#inputs: covx (covariates, n by p matrix); rho (tunning parameter, scale);
#kernel.index ("e": exponential kernel,"g": gaussion Kernel,"l": linear kernel)
gen.ker<-function(covx,rho,kernel.index)
{
  n <- nrow(covx)
  ker <- matrix(0,n,n)
  for (i in 1:n)
    for (j in 1:n)
    {
      x <- covx[i,]
      y <- covx[j,]
      if (kernel.index=="g")
      {
        ker[i,j] <- exp(-sum((x-y)^2)/rho) #gaussian kernel
      }
      if (kernel.index=="e")
      {
        ker[i,j] <- exp((-sum(x^2)-3*sum((x-y)^2)-sum(y^2))/rho) #exponential kernel
      }
      if (kernel.index=="l")
      {
        ker[i,j] <- (sum(x*y)/rho)  #linear kernel
      }

    }
  ker
}
##############calculate U test statistic###############
#inputs: y (quntitative response vector, n by 1); kernel.index (kernel function:"e","l","g")
#: covx (covariants, n by p)
kerUTest<-function(y,covx,ker,alpha=0.05)
{
  n <- length(y)
  p <- ncol(covx)
  h <- diag(n)-matrix(1,n,n)/n
  #ker <- gen.ker(covx,p,kernel.index)
  ker0 <- ker #original kernel with zero diagonals
  diag(ker0) <- rep(0,n)
  J <- matrix(1,n,n)
  ker <- ker-J%*%ker0/(n-1)-ker0%*%J/(n-1)+J%*%ker0%*%J/n/(n-1) #centralized kernel
  ker1 <- ker
  diag(ker1) <- rep(0,n) #centralized kernel with zero diagonals

  A <- h%*%ker1%*%h
  hkhk <- A%*%ker1
  hkh.hkh <- A*A
  hk <- h%*%ker1
  sd <- as.numeric(sqrt(var(y)))
  z <- (y-mean(y))/sd
  m4 <- mean(z^4)
  m6 <- mean(z^6)
  sigma.hat <- var(y)
  e1 <- tr(hk)
  e2 <- tr(hkhk)
  mu1 <- 2*e1
  mu2 <- e1*(1+2/(n-1))
  delta <- m4-3
  var0 <- e1^2*(-2/(n-1))+e2*(2-12/(n-1))
  var <- var0+delta*(-e1^2/n+6*e2/n+tr(hkh.hkh))
  test <- t(y-mean(y))%*%ker1%*%(y-mean(y))/(sigma.hat)
  v1 <- sum(diag(ker))/n
  a <- var/2/(n-1)^2/v1
  g <- v1/a


  ####pvalue#######
  pv.norm <- pnorm(test/sqrt(var),lower.tail = F)
  pv.chisq <- pchisq((test/(n-1)+v1)/a,g,lower.tail = F)
  #if (test>qnorm(0.95,0,sqrt(var))) #testing using asymptotic normal distribution
  #re.norm <- 1
  #if ((test/(n-1)+v1)>a*qchisq(0.95,g)) #testing using scaled chiseq distribution
  #re.u <- 1
  ####power estimation of alpha level test#####
  alpha1 <- y-mean(y)
  sigma.h1 <- var(y)
  shift <- t(alpha1)%*%ker1%*%alpha1/sigma.h1
  shift <- shift/sqrt(var)
  power.chisq <- pchisq(qchisq(1-alpha,g)-shift*sqrt(var/(n-1)^2)/a,g,lower.tail=F)
  return(list(pv.norm=pv.norm, pv.chisq=pv.chisq, power.chisq=power.chisq, shift=shift))
}

##function to select the best kernel function from "e", "l" and "g" kernels.
#value returns to a list of location shift, where the best kernel has the largest shift
KerSel <- function(y,covx)
{
  p <- ncol(covx)
  ker.e <- gen.ker(covx,p,"e")
  ker.l <- gen.ker(covx,p,"l")
  ker.g <- gen.ker(covx,p,"g")
  test.e <- kerUTest(y,covx,ker.e)
  shift.e <- test.e$shift
  test.l <- kerUTest(y,covx,ker.l)
  shift.l <- test.l$shift
  test.g <- kerUTest(y,covx,ker.g)
  shift.g <- test.g$shift
  print("the best kernel has the largest shift")
  return(list(shift.e=shift.e,shift.l=shift.l,shift.g=shift.g))
}

