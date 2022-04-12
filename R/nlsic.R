# Copyright 2011-2022, INRAE/INSA/CNRS
# Author: Serguei SOKOL, sokol@insa-toulouse.fr
# License: GPL v2 http://www.gnu.org/licenses/gpl-2.0.html

#' Join elements into a string
#' @description
#' convert elements of vector v (and all following arguments)
#' in strings and join them using sep as separator.
#'
#' @param sep A string used as a separator
#' @param v A string vector to be joined
#' @param ... other variables to be converted to strings and joined
#' @return A joined string
#' @examples
#' join(" ", c("Hello", "World"))
#' @export
join=function(sep, v, ...) {
  return(paste(c(v,...), sep="", collapse=sep))
}

#' Non Linear Least Squares with Inequality  Constraints
#'
#' @description
#' Solve non linear least squares problem  \code{min_par ||r(par,...)$res||}
#' with optional inequality constraints \code{u\%*\%par >= co}
#' and optional equality constraints \code{e\%*\%par = eco}
#' @param par initial values for parameter vector. It can be in non feasible domain,
#' i.e. for which \code{u\%*\%par >= co} is not always respected.
#' @param r a function calculating residual vector
#'   by a call to \code{r(par, cjac, ...)}
#'   where par is a current parameter vector,
#'   and cjac is set to TRUE if Jacobian is required or FALSE if not.
#'   A call to \code{r()} must return a list having a field "res" containing the
#'   residual vector and optionally a field "jacobian"
#'   when cjac is set to TRUE.
#' @param ... parameters passed through to r()
#' @param u constraint matrix in \code{u\%*\%par >= co}
#' @param co constraint right hand side vector
#' @param control a list of control parameters ('=' indicates default values):
#'   - errx=1.e-7 error on l2 norm of the iteration step sqrt(pt*p).
#'   - maxit=100 maximum of newton iterations
#'   - btstart=1 staring value for backtracking coefficient
#'   - btfrac=0.5 (0;1) by this value we diminish the step till tight up
#'     to the quadratic model of norm reduction in backtrack (bt) iterations
#'   - btdesc=0.1 (0;1) how good we have to tight up to the quadratic model
#'     0 - we are very relaxe, 1 - we are very tight (may need many bt iterations)
#'   - btmaxit=15 maximum of backtrack iterations
#'   - btkmin=1.e-7 a floor value for backtracking fractioning
#'   - trace=0 no information is printed during iterations (1 - print is active)
#'   - rcond=1.e10 maximal condition number in QR decomposition
#'   starting from which a matrix is considered as numerically rank
#'   deficient. The inverse of this number is also used as a measure of very
#'   small number.
#'   - ci = list of options relative to confidence interval estimation
#'      + p=0.95 p-value of confidence interval. If you need a plain sd value,
#'        set p-value to 0.6826895
#'      + report=T report (or not and hence calculate or not) confidence
#'       interval information (in the field hci, as 'half confidence interval' width)
#'   - history=TRUE report or not the convergence history
#'   - adaptbt=FALSE use or not adaptive backtracking
#'   - mnorm=NULL a norm matrix for a case sln==TRUE (cf next parameter)
#'   - sln=FALSE use or not (default) least norm of the solution
#'   (instead of least norm of the increment)
#'   - maxstep=NULL a maximal norm for increment step (if not NULL), must be positive
#'   - monotone=FALSE should or not the cost decrease be monotone. If TRUE, then at
#'    first non decrease of the cost, the iterations are stopped with a warning message.
#' @param e linear equality constraint matrix in \code{e\%*\%par = eco} (dense)
#' @param eco right hand side vector in equality constraints
#' @param flsi function solving linear least squares problem. For a custom function,
#' see interfaces in \code{lsi} or \code{lsi_ln} help pages.
#' @return a list with following components (some components can be absent depending
#' on 'control' parameter)
#'    - 'par' estimated values of par
#'    - 'lastp' the last LSI solution during non linear iterations
#'    - 'hci' vector of half-width confidence intervals for par
#'    - 'ci_p' p-value for which CI was calculated
#'    - 'ci_fdeg' freedom degree used for CI calculation
#'    - 'sd_res' standard deviation of residuals
#'    - 'covpar' covariance matrix for par
#'    - 'laststep' the last increment after possible back-tracking iterations
#'    - 'normp' the norm of lastp
#'    - 'res' the last residual vector
#'    - 'prevres' residual vector from previous non linear iteration
#'    - 'jacobian' the last used jacobian
#'    - 'retres' last returned result of \code{r()} call
#'    - 'it' non linear iteration number where solution was obtained
#'    - 'btit' back-tracking iteration number done during the last non linear iteration
#'    - 'history' list with convergence history information
#'    - 'error' error code: 0 - normal end, 1 - some error occurred, see message in 'mes'
#'    - 'mes' textual message explaining what problem was in case of error
#' @details
#' Solving method consist in sequential LSI problems globalized by backtracking technique.
#' If e, eco are not NULL, reduce jacobian to basis of e's kernel before \code{lsi()} call.\cr
#' NB. If the function \code{r()} returns a list having a field "jacobian" it is
#' supposed to be equal to the jacobian dr/dpar.
#' If not, numerical derivation numDeriv::jacobian()
#' is automatically used for its calculation.\cr
#' NB2. nlsic() does not call stop() on possible errors. Instead, 'error' field is set to 1
#' in the returned result. This is done to allow a user to examine the current state
#' of data and possibly take another path then to merely stop the program. So, a
#' user must allways check this field at return from nlsic().\cr
#' NB3. User should test that field 'mes' is not NULL even when error is 0. It may
#' contain a warning message.
#' @seealso [lsi], [lsi_ln], [uplo2uco]
#' @examples
#' # solve min_{a,b} ||exp(a*x+b)-meas||, a,b>=1
#' a_true=1; b_true=2; x=0:5
#' # simulation function
#' sim=function(par, x) exp(par[["a"]]*x+par[["b"]])
#' # residual function to be passed to nlsic()
#' resid=function(par, cjac, ...) {
#'   dots=list(...)
#'   s=sim(par, dots$x)
#'   result=list(res=s-dots$meas)
#'   if (cjac) {
#'      result$jacobian=cbind(a=s*dots$x, b=s)
#'   }
#'   result
#' }
#' # simulated noised measurements for true parameters
#' set.seed(7) # for reproducible results
#' meas=sim(c(a=a_true, b=b_true), x)+rnorm(x)
#' # starting values for par
#' par=c(a=0, b=0)
#' # prepare constraints
#' uco=uplo2uco(par, lower=c(a=1, b=1))
#' # main call: solve the problem
#' fit=nlsic(par, resid, uco$u, uco$co, control=list(trace=1), x=x, meas=meas)
#' if (fit$error == 1) {
#'    stop(fit$mes)
#' } else {
#'    print(fit$par) # a=1.001590, b=1.991194
#'    if (! is.null(fit$mes)) {
#'       warning(fit$mes)
#'    }
#' }
#' @export
nlsic=function(par, r, u=NULL, co=NULL, control=list(), e=NULL, eco=NULL, flsi=lsi, ...) {
#print(match.call())
   n=length(par)
   m=NROW(co)
   nm_par=names(par)

   # predefined controls, then overwritten by actual values
   con=list(errx=1.e-7, maxit=100, btstart=1., btfrac=0.5, btdesc=0.1,
      btmaxit=15, btkmin=1.e-7, trace=0, rcond=1.e10, ci=list(p=0.95, report=T),
      history=FALSE, adaptbt=FALSE, mnorm=NULL, sln=FALSE, maxstep=NULL, monotone=FALSE)
   nmsC=names(con)
   nmsci=names(con$ci)
   # predefined confidence interval parameters, then overwritten
   ci=con$ci
   if (!is.null(control$ci)) {
      ci[(namci <- names(control$ci))] <- control$ci
   } else {
      namci <- names(con$ci)
   }
   con[(namc <- names(control))] <- control
   con$ci=ci
   tol=1.e-10
   mes=NULL
   if (length(noNms <- namc[!namc %in% nmsC])) {
      mes=join("", "nlsic: unknown name(s) in control: ", join(", ", noNms))
      return(list(par=NULL, error=1, mes=mes))
   }
   if (length(noNms <- namci[!namci %in% nmsci])) {
      mes=join("", "nlsic: unknown names in control$ci: ", join(", ", noNms))
      return(list(par=NULL, error=1, mes=mes))
   }
   if (!is.null(con$maxstep) && con$maxstep <= 0.) {
      mes=sprintf("nlsic: maxstep parameter in control list must be positive, got %g instead", con$maxstep)
      return(list(par=NULL, error=1, mes=mes))
   }
   if (!is.null(con$maxstep) && con$maxstep <= con$errx) {
      mes=sprintf("nlsic: maxstep parameter in control list is less than errx (maxstep=%g, errx=%g)", con$maxstep, con$errx)
      return(list(par=NULL, error=1, mes=mes))
   }

   # history matrices
   hpar=hp=hstep=hres=c()
   if (con$history) {
      hpar=cbind(received=par)
      hist=list()
   } else {
      hist=NULL
   }

   if (!is.null(e)) {
      e=matrix(e, ncol=n)
      ne=nrow(e)
   } else {
      ne=0
   }
   econstr=!is.null(e) && !is.null(eco) && nrow(e)>0 && ncol(e)==length(par)

   # make u matrix if u=NULL
   if (is.null(u)) {
      u=matrix(0, nrow=0, ncol=n)
      co=c()
   }
   if (econstr) {
      nr=nrow(e)
      nc=ncol(e)
      if (nr > nc) {
         return(list(par=NULL, error=1, mes="nlsic: matrix e is over determined."))
      }
      # affine transform epar -> par s.t. e*par=eco
      # par=parp+Null(t(e))%*%epar, where epar is new parameter vector
      nte=Nulla(t(e))
      qe=attr(nte, "qr")
      if (qe$rank < nr) {
         return(list(par=NULL, error=1, mes="nlsic: matrix e is rank deficient."))
      }
      # particular solution
      parp=qr.qy(qe, c(backsolve(qe$qr, eco[qe$pivot], transpose=T), double(n-nr)))
      ue=u%*%nte
      epar=crossprod(nte, c(par)-parp)
      par[]=nte%*%epar+parp
   }
   ine=if (nrow(u) > 0) co-u%*%c(par) else NULL
   if (!is.null(ine) && any(ine > 1./con$rcond)) {
      # project par into feasible domain
      # in case if residual function is defined only on feasible domain
      if (econstr) {
         x=ldp(ue, ine, con$rcond)
         if (is.null(x)) {
            return(list(
               par=NULL,
               error=1,
               mes="nlsic: Unfeasible constraints at starting point"
            ))
         }
         laststep=nte%*%x
         par[]=c(par)+laststep
      } else {
         x=ldp(u, ine, con$rcond)
         if (is.null(x)) {
            mes="nlsic: Unfeasible constraints at starting point (u%*%par-co>0):"
            if (!is.null(rownames(u))) {
               mes=join("\n", c(mes, paste(rownames(u), " (", -ine, ")", sep="")))
            }
            return(list(
               par=NULL,
               error=1,
               mes=mes
            ))
         }
         laststep=x
         par[]=c(par)+laststep
      }
   }

   # newton globalized iterations
   it=0
   btit=0
   converged=F
   laststep=par
   laststep[]=0.
   p=laststep
   a=res=lres=vres=normp=NULL
   while (!converged && it < con$maxit) {
      mes=NULL
      if (it == 0) {
#browser()
         # residual vector
#cat("jac it=0", par, sep="\n")
         lres=r(par, cjac=TRUE, ...)
         if (!is.null(lres$err) && lres$err) {
            return(list(
               par=c(par),
               retres=lres,
               it=it,
               btit=btit,
               error=1,
               mes=paste("nlsic: Problem in first residual calculation",
                  lres$mes, sep="\n"))
            )
         }
         res=as.numeric(lres$res)
         iva=which(!is.na(res)) # valid entries in res
         if (length(iva)==0) {
            # we check it only at it=0 hoping that NA-structure remains the same
            # all the long ietrations
            return(list(
               par=NULL,
               error=1,
               mes="nlsic: no valid residual value."
            ))
         }
         vres=res[iva]
         b=-res
         norb=sqrt(norm2(vres))
         #cat("norb0=", norb, "\n", sep="")
         norres=norb
         klast=0
         if (con$trace) {
            cat("it=", it, "\tres=", norb, "\n", sep="")
         }
         if (norb < 1.e-20) {
            # we are already at solution point
            converged=TRUE
            normp=0
            a=NULL
            next
         }
      } else {
         b=-res
         norb=norres
      }
      # jacobian
      a=NULL
      if ("jacobian" %in% names(lres) && !is.null(lres$jacobian)) {
         a=lres$jacobian[iva,,drop=FALSE]
         colnames(a)=nm_par
      } else {
         # may be we are here because the last call was with cjac=FALSE
         if (it > 0) {
#cat("jac", par, sep="\n")
            lres=r(par, cjac=TRUE, ...)
            res=as.numeric(lres$res)
            iva=!is.na(res)
            vres=res[iva]
            b=-res
            norres=sqrt(norm2(vres))
            norb=norres
            if ("jacobian" %in% names(lres) && !is.null(lres$jacobian)) {
               a=lres$jacobian[iva,,drop=FALSE]
               colnames(a)=nm_par
            }
         }
      }
      if (is.null(a)) {
         if (requireNamespace("numDeriv", quietly = TRUE)) {
            a=numDeriv::jacobian(function(x) as.numeric(r(x, cjac=FALSE, ...)$res), par)[iva,,drop=FALSE]
            colnames(a)=nm_par
         } else {
            return(list(
               par=c(par),
               it=it,
               btit=btit,
               error=1,
               history=hist,
               jacobian=a,
               retres=lres,
               mes=paste("nlsic: Not found 'jacobian' in results of 'r()' when cjac=TRUE and numDeriv package is not available.",
                       attr(p, "mes"), sep="\n")))
         }
      }
      # solve linear ls
      #browser()
      ine=if (nrow(u) > 0) co-u%*%c(par) else NULL
      # newton direction
      if (econstr) {
         if (con$sln) {
            p=flsi(a%*%nte, -vres, ue, ine, rcond=con$rcond, con$mnorm, -c(par))
         } else {
            p=flsi(a%*%nte, -vres, ue, ine, rcond=con$rcond)
         }
         p=nte%*%p
      } else {
         if (con$sln) {
            p=flsi(a, -vres, u, ine, rcond=con$rcond, con$mnorm, -c(par))
         } else {
            p=flsi(a, -vres, u, ine, rcond=con$rcond)
         }
      }
      stopifnot(all(u%*%(par+p)-co >= -tol))
      names(p)=nm_par
#browser()
      if (con$history) {
         hpar=cbind(hpar, par)
         colnames(hpar)[ncol(hpar)]=it
         hp=cbind(hp, c(p))
         hres=cbind(hres, res)
         hist=list(par=hpar, p=hp, res=hres, step=hstep)
      }
      if (anyNA(p)) {
         names(par)=nm_par
         return(list(
            par=c(par),
            it=it,
            btit=btit,
            error=1,
            history=hist,
            jacobian=a,
            retres=lres,
            mes=paste("nlsic: Problem in solving linear problem with constraint",
               attr(p, "mes"), sep="\n")))
      } else if (! is.null(attr(p, "mes"))) {
         mes=join("\n", mes, attr(p, "mes"))
      }

      normp=sqrt(norm2(p))
      converged=(normp <= con$errx)
      if (normp == 0) {
         ############### no need for backtracking at this stage
         laststep=p
         par[]=c(par)+laststep
         it=it+1
         btit=0
         resdecr=NULL
#cat("normp=0", par, sep="\n")
         res=as.numeric(r(par, cjac=FALSE, ...)$res)
         iva=!is.na(res)
         vres=res[iva]
         norres=sqrt(norm2(vres))
         if (con$history) {
            hstep=cbind(hstep, c(laststep))
            hist=list(par=hpar, p=hp, res=hres, step=hstep)
         }
         if (con$trace && (it < 10 || !it%%10)) {
            cat("it=", it, "\tres=", norres, "\tnormstep=", normp, "\tbtk=", klast, "\n", sep="")
         }
         break
      }
      # check if the direction p is descending
      ap=as.double(a%*%p)
      n2ap=norm2(ap)
      resap=crossprod(vres, ap)
      if (resap > 0.) {
         laststep=p*NA
         mes=paste("nlsic: LSI returned not descending direction",
            attr(p, "mes"), sep="\n")
         break
      }
      k=con$btstart; # fraction of p
      if (!is.null(con$maxstep)) {
         kms=con$maxstep/normp
         if (kms < k) {
            k=kms
         }
      }
      # backtrack iterations
      btit=0
      descending=FALSE
      while (!descending && btit < con$btmaxit && k >= con$btkmin) {
         laststep=k*p
#cat("btit=", btit, par+c(laststep), sep="\n")
         lres=try(r(par+c(laststep), cjac=FALSE, ...))
         if (inherits(lres, "try-error")) {
            k=max(k*con$btfrac, con$btkmin)
            btit=btit+1
            next
         }
         if (isTRUE(lres$err != 0)) {
            return(list(
               par=c(par),
               laststep=c(laststep),
               retres=lres,
               it=it,
               btit=btit,
               error=1,
               mes=paste("nlsic: Problem in residual calculation",
                  lres$mes, sep="\n"))
            )
         }
         res=as.numeric(lres$res)
         iva=!is.na(res) # valid entries in res
         vres=res[iva]
         norres=sqrt(norm2(vres))
         #cat("norres=", norres, "\n", sep="")
         #scaprod=crossprod(vres, ap)-resap
         #resdecr=sqrt(scaprod); # residual decreasing
         realdecr=crossprod(b[iva]-vres, vres+b[iva]) # normally it is positive
         lindecr=-k*(2*resap+k*n2ap)
         #descending=(
         #   scaprod >=0. &&
         #   sqrt(scaprod) >= con$btdesc*sqrt(n2ap*k) &&
         #   norres < norb
         #)
         descending=(realdecr>=con$btdesc*lindecr)
         klast=k
         if (con$adaptbt) {
            h=max(min(2., 1./(2.-(norres-norb)*(norres+norb)/(resap*k))), con$btfrac)
            # check that p*k*h won't violate inequalities (for h > 1 only)
            if (h > 1. && nrow(u) > 0) {
               if (any(u%*%(c(par)+(k*h)*p)-co <= -tol)) {
                  h=con$btfrac
               }
            }
            k=k*h
            #cat("h=", h, "klast=", klast, "k=", k, "\n")
         } else {
            k=k*con$btfrac
         }
         k=max(k, con$btkmin)
         btit=btit+1
      }
      if (inherits(lres, "try-error")) {
         return(list(
            par=c(par),
            laststep=rep(NA, length(par)),
            retres=lres,
            it=it,
            btit=btit,
            error=1,
            mes="nlsic: residual calculation failed during all backtracking iterations"
         ))
      }

      par[]=c(par)+laststep
      it=it+1
      if (con$monotone) {
         if (norres > norb) {
            converged=T
            mes=join("\n", mes, "nlsic: Non monotone descrease of cost function has occured.")
         }
      }
      if (con$trace && ((it < 10 || !it%%10) || converged)) {
         cat("it=", it, "\tres=", norres, "\tnormstep=", normp, "\tbtk=", klast, "\n", sep="")
      }
      if (con$history) {
         hstep=cbind(hstep, laststep)
         hist=list(par=hpar, p=hp, res=hres, step=hstep)
      }
   }
   if (con$trace && !(it < 10 || !it%%10) && !converged) {
      cat("it=", it, "\tres=", norres, "\tnormstep=", normp, "\tbtk=", klast, "\n", sep="")
   }
   if (it >= con$maxit) {
      mes=join("\n", mes, "nlsic: Maximal non linear iteration number is achieved")
   }
   if (btit >= con$btmaxit) {
      mes=join("\n", mes, "nlsic: Maximal backtrack iteration number is achieved")
   }
   # restore names
   names(par)=names(laststep)=names(p)=nm_par
   # confidence interval
   #browser()
   nr=length(vres)
   fdeg=max(1., nr-n+ne) # freedom degrees
   if (con$ci$report) {
      sd_res=sqrt(sum(vres**2)/fdeg)
      if (is.null(a)) {
#cat("first jac", par, sep="\n")
         a=r(par, TRUE, ...)$jacobian
         if (is.null(a)) {
#cat("num jac", par, sep="\n")
            if (requireNamespace("numDeriv", quietly = TRUE)) {
               a=numDeriv::jacobian(function(x) as.numeric(r(x, cjac=FALSE, ...)$res), par)
               colnames(a)=nm_par
            } else {
               return(list(
                  par=c(par),
                  it=it,
                  btit=btit,
                  error=1,
                  history=hist,
                  jacobian=a,
                  retres=lres,
                  mes=paste("nlsic: Not found 'jacobian' in results of 'r()' when cjac=TRUE and numDeriv package is not available.",
                            attr(p, "mes"), sep="\n")))
            }
         }
      }
      if (econstr) {
         jac=a%*%nte
      } else {
         jac=a
      }
      covpar=try(solve(crossprod(jac)), silent=T)
      if (inherits(covpar, "try-error")) {
         sv=svd(jac)
         covpar=tcrossprod(sv$v*rep(1./pmax(sv$d, 1./con$rcond), rep(nrow(sv$v), ncol(sv$v))))
      }
      if (econstr) {
         covpar=nte%*%tcrossprod(covpar, nte)
      }
      dimnames(covpar)=list(nm_par, nm_par)
      d=sqrt(diag(covpar))
      hci=(-stats::qt((1.-con$ci$p)/2., fdeg)*sd_res)*d
      names(hci)=nm_par
   } else {
      sd_res=hci=covpar=ci_fdeg=NULL
   }
   return(list(
      par=c(par),
      lastp=c(p),
      hci=hci,
      ci_p=con$ci$p,
      ci_fdeg=fdeg,
      sd_res=sd_res,
      covpar=covpar,
      laststep=c(laststep),
      normp=normp,
      res=c(res),
      prevres=c(vres),
      jacobian=a,
      retres=lres,
      it=it,
      btit=btit,
      history=hist,
      error=0,
      mes=mes)
   )
}

#' Linear Least Squares with Inequality constraints (LSI)
#'
#' @description
#' solve linear least square problem (min ||A%*%x-b||)
#' with inequality constraints \code{u%*%x>=co}
#' @param a dense matrix A or its QR decomposition
#' @param b right hand side vector. Rows containing NA are dropped.
#' @param u dense matrix of inequality constraints
#' @param co right hand side vector of inequality constraints
#' @param rcond maximal condition number for determining rank deficient matrix
#' @param mnorm dummy parameter
#' @param x0 dummy parameter
#' @return solution vector whose attribute 'mes' may contain a message about possible numerical problems
#' @seealso [lsi_ln], [ldp], [base::qr]
#' @details
#' Method:\cr
#' 1. reduce the problem to ldp (min(xat*xa) => least distance programming)\cr
#' 2. solve ldp\cr
#' 3. change back to x\cr
#' If b is all NA, then a vector of NA is returned.
#'
#' mnrom, and x0 are dummy parameters which are here to make lsi()
#' compatible with lsi_ln() argument list
#' @export
lsi=function(a, b, u=NULL, co=NULL, rcond=1.e10, mnorm=NULL, x0=NULL) {
   tol=1.e-10
   # remove NA from b
   i0=which(is.na(b))
   if (length(i0)) {
      if (length(i0) == length(b)) # all NA
         return(rep(NA_real_, if (is.qr(a)) ncol(a$qr) else ncol(a)))
      if (is.qr(a))
         a=qr.X(a)
      a=a[-i0,,drop=FALSE]
      b=b[-i0]
   }
   if (! is.qr(a)) {
      n=ncol(a)
      aq=base::qr(a, LAPACK=TRUE)
   } else {
      n=ncol(a$qr)
      aq=a
   }
   d=abs(diag(aq$qr))
#cat("d=", d, "\n")
   aq$rank=sum(d>d[1]/rcond)
   if (is.na(aq$rank)) {
      x=rep(NA, n)
      attr(x, "mes")="lsi: Rank could not be estimated in least squares\n"
      return(x)
   }
   if (aq$rank < n) {
      x=rep(NA, n)
      attr(x, "mes")=
         paste("lsi: Rank deficient matrix in least squares\n",
         n-aq$rank,
         " unsolvable variable(s):\n",
         paste(dimnames(a)[[2]][aq$pivot[-(1:aq$rank)]],
         aq$pivot[-(1:aq$rank)], sep="\t", collapse="\n"), "\n",
         sep="")
      return(x)
   }
   #x0=qr.solve(aq, b)
   x0=qr.coef(aq, b)
   if (!is.null(u) && nrow(u)>0) {
      # we do have inequalities
      # prepare variable change
      ut=t(base::backsolve(aq$qr, t(as.matrix(u[, aq$pivot, drop=FALSE])), transpose=TRUE))
      xa=ldp(ut, co-u%*%x0, rcond)
      if (is.null(xa)) {
         # Infeasible constraints detected by ldp
         xa=rep(NA, n)
         attr(xa, "mes")="lsi: ldp() revealed unfeasible constrants"
         return(xa)
      }
      x=xa
      xa=base::backsolve(aq$qr, xa)
      x[aq$pivot]=xa
      x=x+x0
      cou=co-u%*%x
      if (any(cou > tol)) {
         # round errors made fail inequality constraints
         # force them as much as possible even by degrading the residual
         dx=ldp(u, cou, rcond)
         if (!is.null(dx)) {
            x=x+dx
         } # else leave x as is
      }
      #stopifnot(all(u%*%x-co >= -tol))
   } else {
      # plain QR
      x=x0
   }
   return(x)
}

#' Linear Least Squares with Inequality constraints, least norm solution
#'
#' @description
#' solve linear least square problem \code{min_x ||A*x-b||}
#' with inequality constraints \code{u\%*\%x >= co}
#' If A is rank deficient, least norm solution \code{||mnorm\%*\%(x-x0)||} is used.
#' If the parameter mnorm is NULL, it is treated as an identity matrix.
#' If the vector x0 is NULL, it is treated as 0 vector.
#' @param a dense matrix A or its QR decomposition
#' @param b right hand side vector
#' @param u dense matrix of inequality constraints
#' @param co right hand side vector of inequality constraints
#' @param rcond maximal condition number for determining rank deficient matrix
#' @param mnorm norm matrix (can be dense or sparse) for which \code{\%*\%} operation with a dense vector is defined
#' @param x0 optional vector from which a least norm distance is searched for
#' @return solution vector whose attribute 'mes' may contain a message about possible numerical problems
#' @seealso [lsi], [ldp], [base::qr]
#' @export
lsi_ln=function(a, b, u=NULL, co=NULL, rcond=1.e10, mnorm=NULL, x0=NULL) {
#if(anyNA(u))
#   browser()
   mes=""
   tol=1.e-10
   if (! is.qr(a)) {
      a=as.matrix(a)
      n=ncol(a)
      aq=base::qr(a, LAPACK=TRUE)
   } else {
      n=ncol(a$qr)
      aq=a
   }
   d=diag(aq$qr)
   aq$rank=sum(abs(d)>abs(d[1])/rcond)
   if (aq$rank == 0) {
      # matrix a is zero => no x can diminish the residual norm
      if (!is.null(u) && nrow(u) > 0) {
         # there are some inequalities to deal with
         if (any(co > tol)) {
            # at least we satisfy inequalities
            x=ldp(u, co, rcond)
            if (is.null(x)) {
               x=rep(NA, n)
               attr(x, "mes")="lsi_ln: matrix a is zero rank and ldp() revealed unfeasible constrants"
            }
            #stopifnot(all(u%*%x-co >= -tol))
            return(x)
         } else {
            # they are already satisfied => return 0 vector
            return(double(n))
         }
      } else {
         # no inequalities
         return(double(n))
      }
   }
   # here after the rank is > 0
   rdefic=aq$rank < n
   if (!rdefic) {
      # just passthrough params to lsi()
      return(lsi(aq, b, u, co, rcond))
   }
   # prepare free variable substitution
   ndefic=n-aq$rank
   i1=seq(len=aq$rank)
   ndefic=n-aq$rank
   i2=seq(aq$rank+1, n, len=ndefic)
   r=base::qr.R(aq)[i1,,drop=FALSE]
   r1=r[,i1,drop=FALSE]
   r2=r[,i2,drop=FALSE]
   qrr=base::qr(t(r), LAPACK=T)
   # null space and its complement bases
   ba=base::qr.Q(qrr, complete=T)
   bal=ba[,i1,drop=FALSE]
   ban=ba[,i2,drop=FALSE]
   # least norm without constraints
   x=bal%*%base::backsolve(qrr$qr, qr.qty(aq, b)[qrr$pivot], transpose=TRUE)
   #x0=backsolve(r1, bt) # implicitly, free variables are set to zero
   #x=c(x0, rep(0., ndefic)) # we unpivot just before return
   #browser()
   if (!is.null(u) && nrow(u)>0) {
      up=u[,aq$pivot,drop=FALSE]
      cou=co-up%*%x
      if (any(cou > tol)) {
         uz=up%*%ban
         # minimize ||y;z||
         ut=t(base::backsolve(qrr$qr, t(as.matrix(up%*%bal)), transpose=FALSE))
         ut=cbind(ut, uz)
         ya=ldp(ut, cou, rcond)
         if (!is.null(ya)) {
            # move along the borders to diminish ||a*x-b||
            # the solution is searched for in the form
            # y=ya+bau*yz
            # ia are indexes of border equations to add as equality constraints (active set)
            cou=as.numeric(cou-ut%*%ya)
            ia=which(cou > -tol) # must not be empty if we are here
            if (length(ia) == 0L) {
               # but it can happen for badly conditioned ut
               # let take the inequality closest to 0
               ia=which.max(cou)
            }
            bau=Nulla(t(ut[ia,,drop=FALSE]))
#if (anyNA(ut[-ia,,drop=FALSE]%*%bau))
#   cat("sent NA in lsi_ln(.., u, ...)\n")
            yz=lsi_ln(bau[i1,,drop=FALSE], -ya[i1], ut[-ia,,drop=FALSE]%*%bau, cou[-ia], rcond=rcond)
            if (!anyNA(yz)) {
               # got a solution
               ya=ya+bau%*%yz
               dx=bal%*%base::backsolve(qrr$qr, ya[i1], transpose=T)+ban%*%ya[i2]
               x=x+dx
               attr(x, "mes")=attr(yz, "mes")
            } else {
               x=rep(NA, n)
               attr(x, "mes")=join("\n", attr(x, "mes"), "lsi_ln: failed to glide along borders.")
               return(x)
            }
         } else {
            # unfeasible constraints
            x=rep(NA, n)
            attr(x, "mes")="lsi_ln: unfeasible constraints"
            return(x)
         }
         # solve lsi_ln within just Null space (without modifying ||y||)
         cou=co-up%*%x
         if(is.null(mnorm)) {
            mban=ban
            if (is.null(x0)) {
               x0=-x
            } else {
               x0=x0-x
            }
         } else {
            mban=mnorm%*%ban
            if (is.null(x0)) {
               x0=-mnorm%*%x
            } else {
               x0=mnorm%*%(x0-x)
            }
         }
         z=lsi_ln(mban, x0, uz, cou)
         if (!anyNA(z)) {
            x=x+ban%*%z
            attr(x, "mes")=join("\n", attr(x, "mes"), attr(z, "mes"))
         } else {
            #print(list(a=a, b=b, u=u, co=co))
#browser()
            attr(x, "mes")=join("\n", attr(x, "mes"), "lsi_ln: Unfeasible constraints in null space. (It must not be. Round off errors struck again.)")
         }
      }
   }
   x[aq$pivot]=x
   if (!is.null(u)) {
      cou=co-u%*%x
      if (any(cou > tol)) {
         # round errors made fail inequality constraints
         # force them as much as possible even by degrading the residual
         dx=ldp(u, cou, rcond)
         if (!is.null(dx)) {
            x=x+dx
         } # else leave x as is
      }
      #stopifnot(all(u%*%x-co >= -tol))
   }
   attr(x, "mes")=join("\n", attr(x, "mes"), paste(
      "lsi_ln: Rank deficient matrix in least squares\n",
      ndefic, " free variable(s):\n",
      paste(dimnames(a)[[2]][aq$pivot[i2]],
      aq$pivot[-i1], sep="\t", collapse="\n"),
      "\nLeast L2-norm solution is provided.",
      sep=""))
   return(x)
}

#' Least Distance Problem
#'
#' @description
#' Solve least distance programming: find x satisfying \code{u\%*\%x >= co} and s.t. min(||x||)
#' by passing to nnls (non negative least square) problem.
#' @param u a dense matrix of inequality constraints
#' @param co a right hand side vector of inequality constraints
#' @param rcond maximal condition number for determining rank deficient matrix
#' @return solution vector or NULL if constraints are unfeasible
#' @importFrom nnls nnls
#' @export
ldp=function(u, co, rcond=1.e10) {
   m=NROW(u)
   n=ncol(u)
   tol=1.e-10
   if (m == 0) {
      # no inequality to satisfy => trivial case
      return(double(n))
   }
   rcond=abs(rcond)
   # eliminate NA from co
   i0=which(is.na(co))
   if (length(i0) > 0) {
      u=u[-i0,,drop=FALSE]
      co=co[-i0]
      m=m-length(i0)
      if (m == 0) {
         # no inequality to satisfy => all NA
         return(rep(NA_real_, n))
      }
   }
   # eliminate 0 rows from u
   maxu=apply(u, 1, function(v) max(abs(v)))
   i0=which(maxu < 1.e-13)
   if (length(i0) > 0) {
      if (all(co[i0] < tol)) {
         u=u[-i0,,drop=FALSE]
         co=co[-i0]
         m=m-length(i0)
         if (m == 0) {
            # no inequality to satisfy => trivial case
            return(double(n))
         }
      } else {
         # unfeasible 0 constraints
         return(NULL)
      }
   }

   maxco=max(co)
   if (maxco<=tol) {
      # all rhs are <= 0 => trivial case
      return(double(n))
   }
   e=rbind(t(u), t(co))
   f=double(n+1L)
   f[n+1]=1.
   resnnls=nnls::nnls(e, f)
   feasible=sqrt(resnnls$deviance) > tol && resnnls$residuals[n+1] != 0.
   if (feasible) {
      x=resnnls$residuals[1:n]/(-resnnls$residuals[n+1])
      # check for numerical stability problems
      ux=u%*%x
      cou=co-ux
      if (FALSE && any(cou > tol)) { # FALSE because it can worsen the situation
         # second trial
         e[n+1,]=cou
         rn=nnls::nnls(e, f)
         if (rn$residuals[n+1]!=0.) {
            x=x-rn$residuals[1:n]/rn$residuals[n+1]
         } else {
            x=NULL
         }
      } else if (all(cou<0) && !all(x==0.)) {
         # round off error pushed the solution inside of feasible domain
         # shorten it till the closest border
         alpha=ux/co
         alpha=min(alpha[alpha>=1.], na.rm=T)
         x=x/alpha
      }
      #stopifnot(all(u%*%x-co >= -tol))
   } else {
      x=NULL
   }
   return(x)
}
norm2=function(v)crossprod(as.numeric(v))[1L]
ls_ln_old=function(a, b, rcond=1.e10) {
   # least squares with least norm
   # no LAPACK in the second qr()
   if (! is.qr(a)) {
      a=base::qr(a, LAPACK=TRUE)
   }
   n=ncol(a$qr)
   d=diag(a$qr)
   a$rank=sum(abs(d)>abs(d[1])/rcond)
   rdefic=a$rank < n
   i1=1:a$rank
   i2=(a$rank+1):n
   r=base::qr.R(a)[i1,,drop=FALSE]
   if (!rdefic) {
      # plain ls
      x=base::backsolve(r, qr.qty(a, b))
      x[a$pivot]=x
      return(x)
   }
   # prepare free variable substitution
   r1=r[,i1,drop=FALSE]
   r2=r[,i2,drop=FALSE]
   qrr=qr(t(r), LAPACK=F)
   # null space and its complement bases
   ba=base::qr.Q(qrr, complete=TRUE)
   bal=ba[,i1,drop=FALSE]
   # least norm
   x=bal%*%base::backsolve(qr.R(qrr), qr.qty(a, b), transpose=TRUE)
   x[a$pivot]=x
   return(x)
}

#' Linear Least Squares, least norm solution
#'
#' @param a matrix or its QR decomposition
#' @param b vector of right hand side
#' @param rcond maximal condition number for rank definition
#' @return solution vector
#' @export
ls_ln=function(a, b, rcond=1.e10) {
   # least squares with least norm
   # LAPACK is called in the second qr() too.
   if (! is.qr(a)) {
      a=base::qr(a, LAPACK=TRUE)
   }
   n=ncol(a$qr)
   d=diag(a$qr)
   a$rank=sum(abs(d)>abs(d[1])/rcond)
   if (a$rank == 0) {
      return(double(n)) # return plain 0 vector
   }
   rdefic=a$rank < n
   i1=seq(len=a$rank)
   r=base::qr.R(a)[i1,,drop=FALSE]
   if (!rdefic) {
      # plain ls
      x=base::backsolve(r, qr.qty(a, b))
      x[a$pivot]=x
      return(x)
   }
   ndefic=n-a$rank
   i2=seq(a$rank+1, n, len=ndefic)
   # prepare free variable substitution
   qrr=base::qr(t(r), LAPACK=TRUE)
   # least norm
   x=base::qr.qy(qrr, c(backsolve(qrr$qr, qr.qty(a, b)[qrr$pivot], transpose=TRUE), rep.int(0., ndefic)))
   x[a$pivot]=x
   return(x)
}

#' Linear Least Squares problem with inequality and equality constraints, least norm solution
#'
#' @description
#' solve linear least square problem (min ||A%*%x-b||)
#' with inequality constraints u%*%x>=co and equality constraints e%*%x=ce
#' Method:
#' reduce the pb to lsi_ln on the null-space of e
#' @param a dense matrix A or its QR decomposition
#' @param b right hand side vector
#' @param u dense matrix of inequality constraints
#' @param co right hand side vector of inequality constraints
#' @param e dense matrix of equality constraints
#' @param ce right hand side vector of equality constraints
#' @param rcond maximal condition number for determining rank deficient matrix
#' @return solution vector whose attribute 'mes' may contain a message about possible numerical problems
#' @seealso [lsi_ln]
#' @export
lsie_ln=function(a, b, u=NULL, co=NULL, e=NULL, ce=NULL, rcond=1.e10) {
   if (is.qr(a)) {
      n=ncol(a$qr)
      aqr=TRUE
   } else {
      n=ncol(a)
      aqr=FALSE
   }
   if (!is.null(e) && nrow(e) > 0) {
      # x is searched in a form xp+bn*z
      # where xp is a particular solution of e*x=ce
      # bn is a basis of e null space and z is a
      # new unknown vector.
      nr=nrow(e)
      nc=ncol(e)
      if (nr > nc) {
         x=rep(NA, n)
         attr(x, "mes")="lsie_ln: matrix e is over determined."
         return(x)
      }
      bn=Nulla(t(e))
      qe=attr(bn, "qr")
      if (qe$rank < nr) {
         x=rep(NA, n)
         attr(x, "mes")="lsie_ln: matrix e is rank deficient."
         return(x)
      }
      # particular solution
      xp=base::qr.qy(qe, c(backsolve(qe$qr, ce[qe$pivot], transpose=T), double(n-nr)))
      if (aqr) {
         R=base::qr.R(a, complete=TRUE)
         b=b-base::qr.qy(a, R%*%xp[a$pivot])
         a=base::qr.qy(a, R%*%bn[a$pivot,,drop=FALSE])
      } else {
         b=b-a%*%xp
         a=a%*%bn
      }
      if (!is.null(u)) {
         co=co-u%*%xp
         u=u%*%bn
      }
      z=lsi_ln(a, b, u, co, rcond)
      x=xp+bn%*%z
      if (!is.null(attr(z, "mes"))) {
         attr(x, "mes")=attr(z, "mes")
      }
      x
   } else {
      lsi_ln(a, b, u=u, co=co, rcond=rcond)
   }
}

#' Linear Least Squares, least norm solution (by svd)
#'
#' @description
#' Least squares \code{a\%*\%x ~= b} of least norm ||x|| by using svd(a)
#' @param a dense matrix
#' @param b right hand side vector
#' @param rcond maximal condition number for determining rank deficient matrix
#' @return solution vector
#' @export
ls_ln_svd=function(a, b, rcond=1.e10) {
   sa=svd(a)
   d=sa$d
   rank=sum(d > d[1L]/rcond)
   if (rank == 0L) {
      return(rep(0., ncol(a)))
   }
   i=seq_len(rank)
   x=sa$v[,i]%*%(crossprod(sa$u[,i], b)/d[i])
   return(x)
}

#' Total Least Squares \code{a\%*\%x ~= b}
#'
#' @param a matrix
#' @param b right hand side vector
#' @return solution vector
#' @export
tls=function(a, b) {
   # total least square by svd
   sab=svd(cbind(a, b))
   n=ncol(a)
   v=sab$v[,length(sab$d)]
   return(v[seq_len(n)]/(-v[n+1L]))
}

#' Regularized Linear Least Squares
#'
#' @description
#' solve linear least square problem \code{(min_x ||a*x-b||)}
#' with inequality constraints ux>=co
#' If a is rank deficient, regularization term \code{lambda^2*||mnorm*(x-x0)||^2}
#' is added to \code{||a*x-b||^2}.
#' @details
#' The rank of a is estimated as number of singular values
#' above \code{d[1]*1.e-10} where \code{d[1]} is the highest singular value.
#' The scalar lambda is an positive number
#' and is calculated as \code{d[1]/sqrt(rcond)} ('rcond' parameter is preserved
#' for compatibility with others lsi_...() functions).
#' At return, lambda can be found in attributes of the returned vector x.
#' NB. lambda is set to NA
#'  - if rank(a)==0 or a is of full rank
#'  - or if there is no inequality.
#' If the matrix mnorm is NULL, it is supposed to be an identity matrix.
#' If the vector x0 is NULL, it is treated as 0 vector.
#' @param a dense matrix A or its QR decomposition
#' @param b right hand side vector
#' @param u dense matrix of inequality constraints
#' @param co right hand side vector of inequality constraints
#' @param rcond used for calculating \code{lambda=d[1]/sqrt(rcond)} where \code{d[1]} is maximal singular value of a
#' @param mnorm norm matrix (can be dense or sparse) for which %*% operation with a dense vector is defined
#' @param x0 optional vector from which a least norm distance is searched for
#' @return solution vector whose attribute 'mes' may contain a message about possible
#'  numerical problems and 'lambda' is regularization parameter used in solution.
#' @seealso [lsi_ln]
#' @export
lsi_reg=function(a, b, u=NULL, co=NULL, rcond=1.e10, mnorm=NULL, x0=NULL) {
   mes=""
   lambda=1./sqrt(rcond)
   tol=1.e-10
   if (is.qr(a)) {
      # get back the matrix from qr
      a=base::qr.X(a)
   }
   n=ncol(a)
   m=nrow(a)
   if (m < n) {
      # make nrow(a) >= ncol(a) by repeating a and b as many times as needed
      nrep=ceiling(n/m)
      a=aperm(array(a, dim=c(m, n, nrep)), c(1,3,2))
      dim(a)=c(m*nrep, n)
      b=rep(b, nrep)
   }
   sva=svd(a)
   d=sva$d
   arank=sum(d > d[1]*tol)
   x0=if (is.null(x0)) double(n) else x0
   if (arank == 0) {
      # matrix a is zero => no x can diminish the residual norm
      if (!is.null(u) && nrow(u) > 0) {
         # there are some inequalities to deal with
         co0=co-u%*%x0
         if (any(co0 > tol)) {
            # at least we satisfy inequalities
            if (is.null(mnorm)) {
               dx=ldp(u, co0)
               if (is.null(dx)) {
                  x=rep(NA, n)
                  attr(x, "mes")="lsi_reg: matrix a is zero rank and ldp() revealed unfeasible constrants"
               }
            } else {
               x=lsi_reg(mnorm, mnorm%*%x0, u, co, rcond, NULL, NULL)
            }
         } else {
            # inequalities are already satisfied
            x=x0
         }
      } else {
         # no inequalities
         x=x0
      }
      attr(x, "lambda")=NA
      return(x)
   }
   # hereafter, the arank is > 0
   i=seq_len(arank)
   if (NROW(u) == 0) {
      # no inequality
      # plain ls_ln
      if (arank == n) {
         x=sva$v[,i]%*%(crossprod(sva$u[,i,drop=FALSE], b)/d[i])
         attr(x, "lambda")=NA
      } else {
         mes=sprintf("lsi_reg: Rank deficient matrix in least squares
There is (are) %d free variable(s).
Regularized L2-norm solution is provided.", n-arank)
         if (is.null(mnorm)) {
            mnorm=diag(1., n)
         }
         lambda=d[1]*lambda
         lam2=lambda**2
         d2=diag(d**2, n)
         mv=mnorm%*%sva$v
         mv2=crossprod(mv)
         bx0=crossprod(mv, mnorm%*%x0)
         bt=d*crossprod(sva$u, b)+lam2*x0
         x=sva$v%*%solve(d2+lam2*mv2, d*bt)
         attr(x, "lambda")=lambda
         attr(x, "mes")=mes
      }
   } else {
      # there are inequalities
      if (arank == n) {
         # a is of full rank => no lambda to use, plain lsi
         invd=1./d
         bt=crossprod(sva$u, b)
         ut=(u%*%sva$v)*rep(invd, rep(n, n))
         y=ldp(ut, co-ut%*%bt)
         x=sva$v%*%(invd*(y+bt))
         attr(x, "lambda")=NA
      } else {
         mes=sprintf("lsi_reg: Rank deficient matrix in least squares
There is (are) %d free variable(s).
Regularized L2-norm solution is provided.", n-arank)
         if (is.null(mnorm)) {
            mnorm=diag(1., n)
         }
         # transformed rhs
         mv=mnorm%*%sva$v
         bx0=crossprod(mv, mnorm%*%x0)
         bt=d*crossprod(sva$u, b)
         # semi-transformed u (it lacks a factor 1/(d**2+(lambda*mv)**2) )
         uv=u%*%sva$v
         mv2=crossprod(mv)
         # ldp to find approximate solution
         lambda=d[1]*lambda
         lam2=lambda**2
         d2=diag(d**2, n)+lam2*mv2
         ut=t(solve(t(d2), t(uv)))
         btt=bt+lam2*bx0
         y=ldp(ut, co-ut%*%btt)
         x=sva$v%*%solve(d2, btt+y)
         attr(x, "lambda")=lambda
         attr(x, "mes")=mes
      }
   }
   return(x)
}

#' Transform box-type inequalities into matrix and vector form
#'
#' @description
#' Transform a set of inequalities  \code{param["name"] >= lower["name"]} and
#'  \code{param["name"] <= upper["name"]} into a list with matrix u and vector co such that
#' u%*%param>=co. In addition to box inequalities above, user can provide linear
#' inequalities in a form like \code{"a+2*c+3*b >= 0"} where 'a', 'b' and 'c' must be names of
#' param components. Numeric and symbolic coefficients and right hand sides
#' are allowed in these expressions. However, symbols must be defined at the moment
#' of calling \code{uplo2uco()} so that expressions containing such symbols could
#' be \code{eval()}-ed to numerical values. All inequalities must be written with '>='
#' sign (not with '<=', '>', ...).
#' @param param a named vector whose names will be used for parsing inequalities
#' @param upper a named numeric vector of upper limits
#' @param lower a named numeric vector of lower limits
#' @param linear a string vector of linear inequalities
#' @return a list with numeric matrix 'u' and vector 'co' such that \code{u%*%param-co>=0}
#' @seealso [equa2vecmat] for parsing linear expressions
#' @export
uplo2uco=function(param, upper=NULL, lower=NULL, linear=NULL) {
   nm_par=names(param)
   nlo=length(lower)
   nup=length(upper)
   nli=length(linear)
   u=matrix(0., nrow=nup+nlo+nli, ncol=length(param))
   co=numeric(nrow(u))
   colnames(u)=nm_par
   rownames(u)=c(
      if (nup) paste(names(upper), " <= ", upper, sep="") else NULL,
      if (nlo) paste(names(lower), " >= ", lower, sep="") else NULL,
      linear
   )
   names(co)=rownames(u)
   # fill u and co
   if (nup) {
      u[seq_len(nup), names(upper)]=diag(-1, nup)
      co[seq_len(nup)]=-upper
   }
   u[nup+seq_len(nlo), names(lower)]=diag(1, nlo)
   co[nup+seq_len(nlo)]=lower

   # parse inequalities
   if (nli) {
      vema=try(equa2vecmat(nm_par, linear, sep=">="))
      if (inherits(vema, "try-error")) {
         stop("Error in parsing inequalities")
      }
      i=nup+nlo+seq_len(nli)
      u[i, ]=vema[,-1L]
      co[i]=vema[,1L]
   }
   return(list(u=u, co=co))
}

#' Parse linear equations/inequalities
#'
#' @description
#' parse a text vector of linear equations and produce a corresponding
#' matrix and right hand side vector
#' @param nm_par a string vector of variable names. It will be used in the symbolic
#' derivation.
#' @param linear string vector of linear equations like \code{"a+2*c+3*b = 0"}
#' @param sep separator of two parts of equations. Use for example
#'    ">=" for linear inequalities
#' @return an augmented matrix. Its first column is the rhs vector.
#' Other columns are named by nm_par. If the vector linear is NULL or its content
#' is empty a NULL is returned
#' @examples
#' equa2vecmat(c("a", "b", "c"), "a+2*c+3*b = 0", "=")
#' @export
equa2vecmat=function(nm_par, linear, sep="=") {
   # Sanity check
   stopifnot(length(sep)==1 && nchar(sep) > 0)
   stopifnot(length(nm_par)>=1 && all(nchar(nm_par)) > 0)

   vlin=sapply(linear, function(it) strsplit(it, sep)[[1]])
   if (length(vlin) && !is.matrix(vlin)) {
      stop(sprintf("Linear (in)equalities are expected to have '%s' sign in them", sep))
   }
   if (length(vlin) == 0) {
      # we are set
      return(NULL)
   }
   # each column in vlin is an (in)equality
   # first row is left part, the second row is the right part
   # derive left and right parts to get the matrix
   ze=double(length(nm_par))
   names(ze)=nm_par
   ze=as.list(ze)
   de=apply(vlin, 2L, function(ineq) {
      le=with(ze, eval(stats::deriv(parse(text=ineq[1]), nm_par)))
      ri=with(ze, eval(stats::deriv(parse(text=ineq[2]), nm_par)))
      v=c(ri-le, attr(le, "gradient")-attr(ri, "gradient"))
      return(v)
   })
   rownames(de)=c("rhs", nm_par)
   return(t(de))
}

lsi_lim=function(a, b, u=NULL, co=NULL, rcond=1.e10, mnorm=NULL, x0=NULL) {
   if (requireNamespace("limSolve", quietly = TRUE))
      suppressWarnings(limSolve::lsei(A=a, B=b, G=u, H=co)$X)
}

#' Null-space basis
#'
#' @description
#' use Lapack for null space basis
#' (derived from MASS::Null)
#' @param M matrix such that \code{t(M)\%*\%B=0} where B is a basis of t(M)'s kernel (aka Null-space)
#' @param rcond maximal condition number for rank definition
#' @return numeric matrix whose columns are basis vectors. Its attribute 'qr' contains QR decomposition of M.
#' @examples
#' Nulla(1:3)
#' @seealso [MASS::Null]
#' @export
Nulla=function (M, rcond=1.e10) {
    tmp <- qr(as.matrix(M), LAPACK=TRUE)
    d=abs(diag(tmp$qr))
    n=length(d)
#browser()
    if (d[1L]==0.) {
       tmp$rank = 0
    } else {
       tmp$rank = sum(d/d[1L] > 1./rcond)
    }
    set <- if (tmp$rank == 0L)
        1L:n
    else
        -(1L:tmp$rank)
    bn=qr.Q(tmp, complete = TRUE)[, set, drop = FALSE]
    attr(bn, "qr")=tmp
    return(bn)
}
