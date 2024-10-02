cctest <- function(formula, data=NULL, df=formula[-2L], ..., tol=1e-7) {
  # Define QR decomposition with row reordering and rank computation:
  QR <- function(x,tol,r=0L,o=c(n,n)[r+n]) {n<-seq_len(nrow(x))
    s<-.colSums(x^2,nrow(x),ncol(x)); s[!s]<-1; x<-x*tcrossprod(n>r,1/sqrt(s))
    q<-qr(x[o,,drop=FALSE],LAPACK=TRUE); t<-abs(diag(q$qr))>tol
    q$rank<-sum(t); q$qr<-q$qr[,t,drop=FALSE]; q$qraux<-q$qraux[t]
    q$o<-o; q$d<-qr.qty(q,x[o,q$pivot[!t],drop=FALSE])*(n>q$rank); q}
  Q <- function(q,y) {y[q$o,]<-qr.qy(q,y);rownames(y)[q$o]<-rownames(q$qr); y}

  # Prepare variables x, y and optional argument as matrices:
  vars <- lapply(list(.a0=df[[length(df)]], .x=formula[[2L]][[3L]],
    .y=formula[[2L]][[2L]]), function(f) do.call(model.matrix.lm,
    list(substitute(~0+f,list(f=f)), data, na.action=NULL),,parent.frame(3)))
  j <- rep.int(0:2, do.call(c,lapply(vars,ncol)))

  # Center rotated variables x, y by removing effects of a:
  cl <- match.call(); cl$df <- cl$tol <- NULL
  cl$formula <- substitute(cbind(.a0,.x,.y)~0+a, list(a=formula[[3L]]))
  mf <- {cl[[1L]]<-quote(stats::model.frame); eval(cl,vars,parent.frame())}
  a <- model.matrix(attr(mf,"terms"), mf); xy <- mf[[1L]]; n <- nrow(xy)
  if (!is.null(h<-model.offset(mf))) xy[,j>0] <- xy[,j>0,drop=FALSE]-h
  if (is.null(w<-model.weights(mf))) w <- rep.int(1,n)
  z <- sqrt(w) + (sqr0<-.Machine$double.xmin^.75); stopifnot(sqr0^2==0)
  a <- a*z; xy <- xy*z
  qa <- QR(a,tol,,order(w,decreasing=TRUE)); ra <- qa$rank
  x <- qr.qty(qa,xy[qa$o,j==1,drop=FALSE])
  y <- qr.qty(qa,xy[qa$o,j==2,drop=FALSE])

  # Compute QR decompositions QxRx and QyRy of the centered data matrices:
  qx <- QR(x,tol,ra); rx <- qx$rank; Qx <- Q(qx,diag(,n,rx))
  qy <- QR(y,tol,ra); ry <- qy$rank; Qy <- Q(qy,diag(,n,ry))

  # Determine residual degrees of freedom (weights are numbers of trials):
  r <- sum(w) - QR(xy[,j==0,drop=FALSE],tol,,qa$o)$rank

  # Compute singular value decomposition of Qy*Qy and new rotated variables:
  SVD <- if (rx && ry) svd(crossprod(Qx,Qy), rx, ry) else
    list(d=numeric(), u=diag(rx), v=diag(ry))
  x. <- Q(qx, rbind(sqrt(r)*SVD$u, matrix(0,n-rx,rx)))
  y. <- Q(qy, rbind(sqrt(r)*SVD$v, matrix(0,n-ry,ry)))

  # Check computability for rows with w=0 (optional, with +sqr0 above):
  dfct <- function(d) .rowSums(abs(d)>tol*z,n,ncol(d)) > 0
  zx <- z; zx[dfct(Q(qa,cbind(qa$d,Q(qx,qx$d))))] <- NaN
  zy <- z; zy[dfct(Q(qa,cbind(qa$d,Q(qy,qy$d))))] <- NaN

  # Compute results:
  V <- sum(SVD$d^2)    # Pillai's statistic
  s <- length(SVD$d); t <- rx*ry; a <- attr(mf,"na.action")
  structure(class="htest", list(
    x = naresid(a,Q(qa,x.)/zx),  # new transformed variables
    y = naresid(a,Q(qa,y.)/zy),
    xinv = crossprod(x.,x)/r,    # inverse coordinate transformations
    yinv = crossprod(y.,y)/r,
    estimate = c(cor=SVD$d),     # canonical correlations (non-negative)
    statistic = c(               # approximate p-values
      "p-value (chi\u00b2 approx.)"=pchisq(V*r, t, lower.tail=FALSE),
      `p-value (F approx.)`=pbeta(V/s, t/2, (r*s-t)/2, lower.tail=FALSE)),
    df.residual = r,             # residual degrees of freedom
    method = "cctest",
    data.name = deparse(substitute(formula), nlines=1L)
  ))
}
