cctest <- function(formula, data=NULL, df=formula[-2L], ..., tol=1e-7,
    stats=FALSE) {
  # Define QR decomposition with row reordering and rank computation:
  QR <- function(x,tol,r=0L,o=c(n,n)[r+n]) {n<-seq_len(nrow(x))
    s<-.colSums(x^2,nrow(x),ncol(x)); s[!s]<-1; x<-x*tcrossprod(n>r,1/sqrt(s))
    q<-qr(x[o,,drop=FALSE],LAPACK=TRUE); t<-abs(diag(q$qr))>tol
    q$rank<-sum(t); q$qr<-q$qr[,t,drop=FALSE]; q$qraux<-q$qraux[t]
    q$o<-o; q$d<-qr.qty(q,x[o,q$pivot[!t],drop=FALSE])*(n>q$rank); q}
  Q <- function(q,y) {y[q$o,]<-qr.qy(q,y)
    rownames(y)<-NULL; rownames(y)[q$o]<-rownames(q$qr); y}

  # Define function for converting data to numeric matrices with same nrow:
  matrices <- function(...) {
    l <- lapply(list(...), function(x) {
      if (is.null(x))
        x <- matrix(numeric(0), 1L)
      if (!is.array(x)) {
        colnms <- levels(x)
        if (is.null(colnms) && is.character(x))
          colnms <- sort.int(unique(x), method="radix")
        if (!is.null(colnms)) {
          x <- diag(length(colnms))[match(x,colnms),,drop=FALSE]
          dimnames(x) <- list(NULL, colnms)
        } else {
          x <- matrix(x)
        }
      }
      storage.mode(x) <- "double"
      x
    })
    is.matrix <- vapply(l, is.matrix, NA)
    nrow <- vapply(l, nrow, 0L)
    n <- nrow[nrow!=1L][1L]
    stopifnot(is.matrix, nrow==1L | nrow==n)
    if (!is.na(n)) for (i in seq_along(l)[nrow==1L])
      l[[i]] <- l[[i]][rep.int(1L,n),,drop=FALSE]
    makenms <- vapply(l, function(x) ncol(x)==1L && is.null(dimnames(x)), NA)
    if (any(makenms)) {
      expr <- substitute(list(...))[-1L]
      for (i in seq_along(l)[makenms])
        dimnames(l[[i]]) <- list(NULL, deparse(expr[[i]], nlines=1L))
    }
    l
  }

  # Prepare variables as matrices:
  f <- list(Y=formula[[2L]][[2L]], X=formula[[2L]][[3L]],
    A=formula[[3L]], A0=df[[length(df)]])
  cl <- match.call(); cl$df <- cl$tol <- cl$stats <- NULL
  if (stats) {   # 'stats' formula syntax (using model.frame, model.matrix)
    cl$formula <- formula
    cl$formula[[2L]] <- substitute(Y+X+A+A0, f); cl$formula[[3L]] <- NULL
    mf <- {cl[[1L]]<-quote(stats::model.frame); eval.parent(cl)}
    vars <- lapply(f, function(f) do.call(model.matrix,
      list(substitute(~0+f,list(f=f)), mf),,parent.frame(3)))
    if (!is.null(h<-model.offset(mf))) {vars$X<-vars$X-h; vars$Y<-vars$Y-h}
    if (is.null(w<-model.weights(mf))) w <- rep.int(1,nrow(vars$A))
    naadjust <- function(x) naresid(attr(mf,"na.action"), x)
  } else {       # simplified syntax (using function matrices)
    cl$formula <- cl$data <- NULL
    env <- new.env(parent=environment(formula))
    assign(envir=env, "|", function(...) do.call(cbind,
      c(deparse.level=0, matrices(...))))
    assign(envir=env, ":", function(...) Reduce(function(x, y) {
      nx <- ncol(x); ny <- ncol(y)
      xe <- x[, rep.int(seq_len(nx),rep.int(ny,nx)), drop=FALSE]
      ye <- y[, rep.int(seq_len(ny),nx), drop=FALSE]
      `colnames<-`(xe*ye, paste(sep=":", colnames(xe), colnames(ye)))
    }, matrices(...)))
    vars <- eval(substitute(.(Y=Y,X=X,A=A,A0=A0), c(f,.=matrices)), data, env)
    args <- lapply(cl[-1L], eval, data, parent.frame())
    w <- as.numeric(do.call(complete.cases, vars))
    if (!is.null(wts<-args$weights)) w[] <- wts
    if (!is.null(rownms<-row.names(data)) && length(rownms)==length(w))
      vars <- lapply(vars, `rownames<-`, rownms)
    if (!is.null(ss<-args$subset))
      {vars <- lapply(vars,`[`,ss,,drop=FALSE); w <- w[ss]}
    naadjust <- function(x) x
  }

  # Center rotated variables X, Y by removing effects of A:
  z <- sqrt(w) + (sqrt0<-.Machine$double.xmin^.75); stopifnot(sqrt0^2==0)
  vars <- lapply(vars, `*`, z)
  na <- lapply(vars, function(x) !w & !complete.cases(x))
  for (i in seq_along(vars)) vars[[i]][na[[i]],] <- 0
  zx <- zy <- z; zx[na$A|na$X] <- zy[na$A|na$Y] <- NA
  qa <- QR(vars$A,tol,,order(w,decreasing=TRUE))
  X <- qr.qty(qa,vars$X[qa$o,,drop=FALSE])
  Y <- qr.qty(qa,vars$Y[qa$o,,drop=FALSE])

  # Compute QR decompositions QxRx and QyRy of the centered data matrices:
  n <- length(z)
  qx <- QR(X,tol,qa$rank); k <- qx$rank; Qx <- Q(qx,diag(,n,k))
  qy <- QR(Y,tol,qa$rank); l <- qy$rank; Qy <- Q(qy,diag(,n,l))

  # Compute singular value decomposition of Qx*Qy and new rotated variables:
  SVD <- if (k && l) svd(crossprod(Qx,Qy), k, l) else
    list(d=numeric(), u=diag(k), v=diag(l))
  x <- Q(qx, rbind(SVD$u, matrix(0,n-k,k)))
  y <- Q(qy, rbind(SVD$v, matrix(0,n-l,l)))

  # Check computability for rows with w=0 (optional, with +sqrt0 above):
  dfct <- function(d) .rowSums(abs(d)>tol*z,n,ncol(d)) > 0
  zx[dfct(Q(qa,cbind(qa$d,Q(qx,qx$d))))] <- NaN
  zy[dfct(Q(qa,cbind(qa$d,Q(qy,qy$d))))] <- NaN

  # Determine residual degrees of freedom (weights are numbers of trials):
  r <- sum(w) - QR(vars$A0,tol,,qa$o)$rank; s <- sqrt(r)

  # Compute results:
  d <- c(cor=SVD$d); t <- k*l; u <- c(beta=r*length(d), gamma=Inf)
  structure(class="htest", list(
    x = naadjust(Q(qa,x)*(s/zx)),  # new transformed variables
    y = naadjust(Q(qa,y)*(s/zy)),
    xinv = (xi<-crossprod(x,X))/s, # inverse coordinate transformations
    yinv = (yi<-crossprod(y,Y))/s,
    estimate = {                   # canonical correlations (non-negative)
      if (t==1L) {p<-crossprod(xi,yi)*d
        if (all(p>0)) d<-c(pos=d); if (all(p<0)) d<-c(neg=d)}; d},
    statistic = c(`p-value`=       # approximate p-values
      replace(pf((1/(sum(d^2)*r)-1/u) / (1/t-1/u), u-t, t), !(0<t & t<u), 1)),
    df.residual = r,               # residual degrees of freedom
    method = "cctest",
    data.name = deparse(formula, nlines=1L)
  ))
}
