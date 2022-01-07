set.seed(7) # for reproducible results
test.join <- function() {
   checkEquals(join(" ", c("Hello", "World"), "!"), "Hello World !")
}
test.Nulla <- function() {
   n=5
   v=rnorm(n)
   base=Nulla(v)
   checkEqualsNumeric(t(v)%*%base, rep(0, n-1))
}
test.lsi <- function() {
   nr=5
   nc=3
   a=matrix(rnorm(nr*nc), nrow=nr, ncol=nc)
   xe=rnorm(nc)
   b=a%*%xe
   # plain ls
   checkEqualsNumeric(lsi(a, b), xe)
   # add inequality
   par=seq_len(nc)
   names(par)=paste("x", par, sep="")
   uco=uplo2uco(par, upper=c(x1=xe[1L]-1))
   # inequality must be enforced
   x=lsi(a, b, uco$u, uco$co)
   checkEqualsNumeric(x[1], xe[1]-1)
   # Kuhn-Tucker condition: t(a)%*%r can be decomposed in t(u) for enforced inequalities
   checkEqualsNumeric((t(a)%*%(a%*%x-b))[-1], rep(0, nc-1))
}
test.lsi_ln <- function() {
   nr=5
   nc=3
   a=matrix(rnorm(nr*nc), nrow=nr, ncol=nc)
   a[,nc]=a[,nc-1] # make rank deficient
   xe=rnorm(nc)
   b=a%*%xe
   # plain ls_ln
   x=lsi_ln(a, b)
   # the last two must be equal
   checkEqualsNumeric(x[nc], x[nc-1])
   # add inequality
   par=seq_len(nc)
   names(par)=paste("x", par, sep="")
   uco=uplo2uco(par, upper=c(x1=xe[1L]-1))
   # inequality must be enforced
   x=lsi_ln(a, b, uco$u, uco$co)
   checkEqualsNumeric(x[1], xe[1]-1)
   # the last two must be equal
   checkEqualsNumeric(x[nc], x[nc-1])
   # Kuhn-Tucker condition: t(a)%*%r can be decomposed in t(u) for enforced inequalities
   checkEqualsNumeric((t(a)%*%(a%*%x-b))[-1], rep(0, nc-1))
}
test.ls_ln <- function() {
   nr=5
   nc=3
   a=matrix(rnorm(nr*nc), nrow=nr, ncol=nc)
   a[,nc]=a[,nc-1] # make rank deficient
   xe=rnorm(nc)
   b=a%*%xe
   # plain ls_ln
   x=ls_ln(a, b)
   # the last two must be equal
   checkEqualsNumeric(x[nc], x[nc-1])
}
test.ls_ln_svd <- function() {
   nr=5
   nc=3
   a=matrix(rnorm(nr*nc), nrow=nr, ncol=nc)
   a[,nc]=a[,nc-1] # make rank deficient
   xe=rnorm(nc)
   b=a%*%xe
   # plain ls_ln
   x=ls_ln_svd(a, b)
   # the last two must be equal
   checkEqualsNumeric(x[nc], x[nc-1])
   # must be consistent with QR solution
   checkEqualsNumeric(x, ls_ln(a, b))
}
test.ldp <- function() {
   # prepare inequality x1 >= 1
   nc=3L
   par=seq_len(nc)
   names(par)=paste("x", par, sep="")
   uco=uplo2uco(par, lower=c(x1=1))
   # solution must be c(1, 0, 0)
   x=ldp(uco$u, uco$co)
   checkEqualsNumeric(x, c(1, 0, 0))
   # x1 >= -1
   uco=uplo2uco(par, lower=c(x1=-1))
   # solution must be c(0, 0, 0)
   x=ldp(uco$u, uco$co)
   checkEqualsNumeric(x, c(0, 0, 0))
   # unfeasible inequalities
   uco=uplo2uco(par, upper=c(x1=-1), lower=c(x1=1))
   x=ldp(uco$u, uco$co)
   checkTrue(is.null(x))
}
test.nlsic <- function() {
   # solve min_{a,b} ||exp(a*x+b)-meas||, a,b>=1
   a_true=1; b_true=2; x=0:5
   # simulation function
   sim=function(par, x) exp(par[["a"]]*x+par[["b"]])
   # residual function to be passed to nlsic()
   resid=function(par, cjac, ...) {
      dots=list(...)
      s=sim(par, dots$x)
      result=list(res=s-dots$meas)
      if (cjac) {
         result$jacobian=cbind(a=s*dots$x, b=s)
      }
      result
   }
   # simulated measurements for true parameters
   meas=sim(c(a=a_true, b=b_true), x)
   # starting values for par
   par=c(a=0, b=0)
   # prepare constraints
   uco=uplo2uco(par, lower=c(a=1, b=1))
   # main call: solve the problem
   fit=nlsic(par, resid, uco$u, uco$co, x=x, meas=meas)
   if (fit$error == 1) {
      stop(fit$mes)
   } else {
      checkEqualsNumeric(fit$par, c(a_true, b_true))
      if (! is.null(fit$mes)) {
         warning(fit$mes)
      }
   }
}
