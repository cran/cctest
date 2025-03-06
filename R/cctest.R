cctest <- function(formula, data=NULL, df=formula[-2L], ..., tol=1e-7) {
  # Define QR decomposition with row reordering and rank computation:
  QR <- function(x,tol,r=0L,o=c(n,n)[r+n]) {n<-seq_len(nrow(x))
    s<-.colSums(x^2,nrow(x),ncol(x)); s[!s]<-1; x<-x*tcrossprod(n>r,1/sqrt(s))
    q<-qr(x[o,,drop=FALSE],LAPACK=TRUE); t<-abs(diag(q$qr))>tol
    q$rank<-sum(t); q$qr<-q$qr[,t,drop=FALSE]; q$qraux<-q$qraux[t]
    q$o<-o; q$d<-qr.qty(q,x[o,q$pivot[!t],drop=FALSE])*(n>q$rank); q}
  Q <- function(q,y) {y[q$o,]<-qr.qy(q,y);rownames(y)[q$o]<-rownames(q$qr); y}

  # Prepare variables as matrices:
  f <- list(Y=formula[[2L]][[2L]], X=formula[[2L]][[3L]],
    A=formula[[3L]], A0=df[[length(df)]])
  cl <- match.call(); cl$df <- cl$tol <- NULL; cl$formula <- formula
  cl$formula[[2L]] <- substitute(Y+X+A+A0, f); cl$formula[[3L]] <- NULL
  mf <- {cl[[1L]]<-quote(stats::model.frame); eval.parent(cl)}
  vars <- lapply(f, function(f) do.call(model.matrix,
    list(substitute(~0+f,list(f=f)), mf),,parent.frame(3)))
  n <- nrow(vars$A)
  if (!is.null(h<-model.offset(mf))) {vars$X <- vars$X-h; vars$Y <- vars$Y-h}
  if (is.null(w<-model.weights(mf))) w <- rep.int(1,n)

  # Center rotated variables X, Y by removing effects of A:
  z <- sqrt(w) + (sqr0<-.Machine$double.xmin^.75); stopifnot(sqr0^2==0)
  vars <- lapply(vars, `*`, z)
  qa <- QR(vars$A,tol,,order(w,decreasing=TRUE)); ra <- qa$rank
  X <- qr.qty(qa,vars$X[qa$o,,drop=FALSE])
  Y <- qr.qty(qa,vars$Y[qa$o,,drop=FALSE])

  # Compute QR decompositions QxRx and QyRy of the centered data matrices:
  qx <- QR(X,tol,ra); rx <- qx$rank; Qx <- Q(qx,diag(,n,rx))
  qy <- QR(Y,tol,ra); ry <- qy$rank; Qy <- Q(qy,diag(,n,ry))

  # Determine residual degrees of freedom (weights are numbers of trials):
  r <- sum(w) - QR(vars$A0,tol,,qa$o)$rank

  # Compute singular value decomposition of Qx*Qy and new rotated variables:
  SVD <- if (rx && ry) svd(crossprod(Qx,Qy), rx, ry) else
    list(d=numeric(), u=diag(rx), v=diag(ry))
  x <- Q(qx, rbind(sqrt(r)*SVD$u, matrix(0,n-rx,rx)))
  y <- Q(qy, rbind(sqrt(r)*SVD$v, matrix(0,n-ry,ry)))

  # Check computability for rows with w=0 (optional, with +sqr0 above):
  dfct <- function(d) .rowSums(abs(d)>tol*z,n,ncol(d)) > 0
  zx <- z; zx[dfct(Q(qa,cbind(qa$d,Q(qx,qx$d))))] <- NaN
  zy <- z; zy[dfct(Q(qa,cbind(qa$d,Q(qy,qy$d))))] <- NaN

  # Compute results:
  V <- sum(SVD$d^2)    # Pillai's statistic
  s <- length(SVD$d); t <- rx*ry; a <- attr(mf,"na.action")
  structure(class="htest", list(
    x = naresid(a,Q(qa,x)/zx),   # new transformed variables
    y = naresid(a,Q(qa,y)/zy),
    xinv = crossprod(x,X)/r,     # inverse coordinate transformations
    yinv = crossprod(y,Y)/r,
    estimate = c(cor=SVD$d),     # canonical correlations (non-negative)
    statistic = c(               # approximate p-values
      "p-value (chi\u00b2 approx.)"=pchisq(V*r, t, lower.tail=FALSE),
      `p-value (F approx.)`=pbeta(V/s, t/2, (r*s-t)/2, lower.tail=FALSE)),
    df.residual = r,             # residual degrees of freedom
    method = "cctest",
    data.name = deparse(substitute(formula), nlines=1L)
  ))
}
