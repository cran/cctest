\name{cctest}
\alias{cctest}
\encoding{UTF-8}
\title{Tests of Independence Based on Canonical Correlations}
\author{Robert Schlicht}
\keyword{htest}
\keyword{multivariate}
\description{
  \code{cctest} estimates canonical correlations between two sets of 
  variables, possibly after removing effects of a third set of variables, and 
  performs a classical multivariate test of (conditional) independence based 
  on Pillai’s statistic.
}
\usage{
  cctest(formula, data=NULL, weights=NULL, ..., tol=1e-07, stats=FALSE)
}
\arguments{
  \item{formula}{An object of the form \code{Y ~ X ~ A}, where \code{Y} 
    represents dependent variables, \code{X} represents a second set of 
    dependent variables or explanatory variables not present under the null 
    hypothesis, and \code{A} represents explanatory variables that remain 
    under the null hypothesis. Typically \code{A} includes at least the 
    constant \code{1} to specify a model with intercepts; unlike 
    \code{\link[stats]{lm}}, the function never adds this automatically.

    If \code{stats} is \code{FALSE} (recommended, see Note), \code{A}, 
    \code{X}, \code{Y} are specified in a simplified notation where all 
    symbols and operators except \code{|} and \code{:} have their regular 
    meaning, with \code{|} used for joining terms (instead of \code{+}; see 
    \code{\link[base]{cbind}}) and \code{:} for multiplying terms, valid 
    anywhere in the expression. The latter differs from \code{*} in that each 
    combination of columns in the arguments contributes a column to the 
    result and multiplication by 0 always yields 0, even for missing values. 
    Single-row arguments in \code{~}, \code{|}, \code{:} are expanded to 
    objects with identical rows, and every \code{factor} or character 
    variable is represented by its full set of indicator (dummy) variables. 
    The values of character variables are \link[base:sort]{sorted} by the 
    \code{"radix"} method.
    %For example, if \code{x} and \code{y} are vectors with elements 
    %\code{5,NA,10} and \code{A,A,B}, then \code{(1|x):y} is the matrix 
    %\code{\tabular{rrrrr}{
    %      \tab 1:A\tab 1:B\tab x:A\tab x:B\cr
    %  [1,]\tab   1\tab   0\tab   5\tab   0\cr
    %  [2,]\tab   1\tab   0\tab  NA\tab   0\cr
    %  [3,]\tab   0\tab   1\tab   0\tab  10
    %}}
  }
  \item{data}{An optional list (or data frame) in which the variables in 
    \code{formula}, \code{weights} and \code{...} are looked for (prior to 
    lookup in \code{|}, \code{:} and in the environment in which the 
    \code{\link[base]{~}} preceding \code{A} was invoked).
  }
  \item{weights}{An optional object of the form \code{~ A0} or \code{w ~ A0}, 
    where \code{A0} replaces \code{A} for the degrees of freedom computation; 
    the default means both are the same. The weights vector \code{w} consists 
    of any nonnegative numbers, equal to 1 by default, that specify how many 
    identical observations each row represents. For rows of data with missing 
    values these numbers are internally replaced with 0, and rows with 
    missing weights are removed.
  }
  \item{\dots}{Additional optional arguments in the case that \code{stats} is 
    \code{TRUE}.
  }
  \item{tol}{The tolerance in the QR decomposition for detecting linear 
    dependencies of the matrix columns.
  }
  \item{stats}{A logical value. If \code{TRUE}, the expressions in 
    \code{formula}, \code{data}, \code{weights} and \code{\dots} are passed 
    to \code{model.frame} and \code{model.matrix} for processing; the 
    \link[stats:formula]{operators} and \link[stats:model.matrix]{expansion} 
    rules defined for the \code{model} part of a \code{formula} object here 
    apply to all parts. If \code{FALSE}, \code{formula} and \code{weights} 
    are interpreted according to the simplified notation.
  }
}
\value{
  A list with class \code{htest} containing the following components:
  \item{x, y}{matrices \eqn{\tilde{X}} and \eqn{\tilde{Y}} of new transformed 
    variables with \eqn{K} and \eqn{L} columns, respectively}
  \item{xinv, yinv}{matrices \eqn{U} and \eqn{V} representing the inverse 
    coordinate transformations}
  \item{estimate}{vector of canonical correlations, i.e., the 
    diagonal elements of \eqn{D}, possibly (only if \eqn{K=L=1}) with a name 
    indicating the direction of the correlation}
  \item{statistic}{vector of p-values based on Pillai’s statistic and 
    classical F (beta) and chi-squared (gamma) approximations}
  \item{df.residual}{the number \eqn{r}}
  \item{method}{the name of the function}
  \item{data.name}{a character string representation of \code{formula} 
    (possibly shortened)}
}

\details{
\code{cctest} unifies various classical statistical procedures that involve 
the same underlying computations, including t-tests, tests in univariate and 
multivariate linear models, parametric and nonparametric tests for 
correlation, Kruskal--Wallis tests, common approximate versions of Wilcoxon 
rank-sum and signed rank tests, chi-squared tests of independence, score 
tests of particular hypotheses in generalized linear models, canonical 
correlation analysis and linear discriminant analysis (see Examples).

Specifically, for the matrices with ranks \eqn{K} and \eqn{L} obtained from 
\eqn{X} and \eqn{Y} by subtracting from each column its orthogonal projection 
on the column space of \eqn{A}, the function computes factorizations 
\eqn{\tilde{X}U} and \eqn{\tilde{Y}V} with \eqn{\tilde{X}} and 
\eqn{\tilde{Y}} having \eqn{K} and \eqn{L} columns, respectively, such that 
both \eqn{\tilde{X}^\top \tilde{X}=rI} and \eqn{\tilde{Y}^\top \tilde{Y}=rI}, 
and \eqn{\tilde{X}^\top \tilde{Y}=rD} is zero except for nonnegative elements 
on its main diagonal arranged in decreasing order. The scaling factor 
\eqn{r}, which should be nonzero, is the dimension of the orthogonal 
complement of the column space of \eqn{A_{0}}.

The function realizes this variant of the singular value decomposition by 
first computing preliminary QR factorizations, taking \eqn{r=1} and without 
the diagonal requirement on \eqn{D}, and then, in a second step, modifying 
these based on a standard singular value decomposition of that matrix. The 
main work is done in a rotated coordinate system where the column space of 
\eqn{A} aligns with the coordinate axes. The basic approach and the rank 
detection algorithm are inspired by the implementations in 
\code{\link[stats]{cancor}} and in \code{\link[stats:lm.fit]{lm}}, 
respectively, although \acronym{LAPACK} routines are used.

The main diagonal elements of \eqn{D}, or singular values, are the estimated 
\dfn{canonical correlations} (Hotelling, 1936) of the variables represented 
by \eqn{X} and \eqn{Y} if these follow a linear model 
\eqn{(X\;\;Y)=A(\alpha\;\;\beta)+(\delta\;\;\epsilon)} with known \eqn{A}, 
unknown \eqn{(\alpha\;\;\beta)} and error terms \eqn{(\delta\;\;\epsilon)} 
that have uncorrelated rows with expectation zero and an identical unknown 
covariance matrix. In the most common case, where \code{A} is given as a 
constant \code{1}, these are the sample canonical correlations (i.e., based 
on simple centering) most often presented in the literature for full column 
ranks \eqn{K} and \eqn{L}. They are always decreasing and between 0 and 1.

In the case of the linear model with independent normally distributed rows 
and \eqn{A_{0}=A}, the ranks \eqn{K} and \eqn{L} equal, with probability 1, 
the ranks of the covariance matrices of the rows of \eqn{X} and \eqn{Y}, 
respectively, or \eqn{r}, whichever is smaller. Under the hypothesis of 
independence of \eqn{X} and \eqn{Y}, given those ranks, the joint 
distribution of the \eqn{J} squared singular values, where \eqn{J} is the 
smaller of the two ranks, is then known and in the case \eqn{r\geq K+L} has a 
probability density (Hsu, 1939, Anderson, 2003, Anderson, 2007) given by
\deqn{\rho (t_{1},...,t_{J})\propto \prod _{j=1}^{J}t_{j}^{(\left|K-L
\right|-1)/2}(1-t_{j})^{(r-K-L-1)/2}\prod _{j'>j}(t_{j}-t_{j'}),} 
\eqn{1\geq t_{1}\geq \cdots \geq t_{J}\geq 0}. For \eqn{J=1} this reduces to 
the well-known case of a single beta distributed \eqn{R^{2}} or equivalently 
an F distributed \eqn{R^{2}/(KL) \over (1-R^{2})/(r-KL)}, with the divisors 
in the numerator and denominator representing the degrees of freedom, or 
twice the parameters of the beta distribution.

Pillai’s statistic is the sum of squares of the canonical correlations, which 
equals, even without the diagonal requirement on \eqn{D}, the squared 
Frobenius norm of that matrix (or trace of \eqn{D^\top D}). Replacing the 
distribution of that statistic divided by \eqn{J} (i.e., of the mean of 
squares) with beta or gamma distributions with first or shape parameter 
\eqn{KL/2} and expectation \eqn{KL/(rJ)} leads to the F and chi-squared 
approximations that the p-values returned by \code{cctest} are based on.

The F or beta approximation (Pillai, 1954, p. 99, p. 44) is usually used with 
\eqn{A_{0}=A} and then is exact if \eqn{J=1}. The chi-squared approximation 
represents Rao’s (1948) score test (with a test statistic that is \eqn{r} 
times Pillai’s statistic) in the model obtained after removing (or 
conditioning on) the orthogonal projections on the column space of 
\eqn{A_{0}} provided that is a subset of the column space of \eqn{A};
see Mardia and Kent (1991) for the case with independent identically 
distributed rows. A more comprehensive derivation of computational 
equivalences of several procedures is provided in a separate 
\href{https://wwwpub.zih.tu-dresden.de/~roschlic/Schlicht2026_correlations.pdf}{manuscript}.
}

\note{
The handling of the weights differs from that in \code{\link[stats]{lm}}. 
Also, to facilitate predictions for rows with zero weights (see Examples), 
the square roots of the weights, used internally for scaling the data, are 
always computed as nonzero numbers, even for zero weights, where they are so 
small that their square is still numerically zero and hence without effect on 
the correlation analysis. For rows of data that imply higher ranks in a 
calculation compensating for the weight, the corresponding rows of the 
returned \code{x} or \code{y} are set to \code{NaN}, although this is visible 
only if they have at least one column.

The simplified \code{formula} notation is intended to provide a simpler, more 
consistent behavior than the legacy \pkg{stats} interface based on 
\code{\link[stats]{terms.formula}}, \code{\link[stats]{model.frame}} and 
\code{\link[stats]{model.matrix}}. Inconsistent or hard-to-predict behavior 
can result in \code{model.matrix}, in particular, from the special 
interpretation of common symbols, the identification of variables by deparsed 
expressions, the locale-dependent conversion of character variables to 
factors and the attempts at reducing linear dependencies subject to 
\code{options("contrasts")}. The \code{stats} option may be removed in the 
future.

For the classical rank tests shown in the Examples, \code{data} must be 
passed to \code{cctest} with subsetting and removal of missing values done 
beforehand.

Running \code{cctest} by itself shows the complete implementation (150 lines).
}

\references{
Hotelling, H. (1936). Relations between two sets of variates. 
  \emph{Biometrika} 28, 321–377.
  \doi{10.1093/biomet/28.3-4.321}, \doi{10.2307/2333955}

Hsu, P.L. (1939). On the distribution of roots of certain determinantal 
  equations. \emph{Annals of Eugenics} 9, 250–258.
  \doi{10.1111/j.1469-1809.1939.tb02212.x}

Rao, C.R. (1948). Large sample tests of statistical hypotheses 
  concerning several parameters with applications to problems of estimation. 
  \emph{Mathematical Proceedings of the Cambridge Philosophical Society} 44, 
  50–57.
  \doi{10.1017/S0305004100023987}

Pillai, K.C.S. (1954). \emph{On some distribution problems in multivariate 
  analysis} (Institute of Statistics mimeo series 88). North Carolina State 
  University, Dept. of Statistics.

Mardia, K.V., Kent, J.T. (1991). Rao score tests for goodness of fit and 
  independence. \emph{Biometrika} 78, 355–363.
  \doi{10.1093/biomet/78.2.355}

Anderson, T.W. (2003). \emph{An introduction to multivariate statistical 
  analysis}, 3rd edition, Ch. 12–13. Wiley.

Anderson, T.W. (2007). Multiple discoveries: distribution of roots of 
  determinantal equations. \emph{Journal of Statistical Planning and 
  Inference} 137, 3240–3248.
  \doi{10.1016/j.jspi.2007.03.008}
}

\seealso{
\code{\link[stats]{cancor}}, \code{\link[stats]{anova.mlm}} in package 
\pkg{stats}; functions in package \pkg{coin}, specifically 
\code{independence_test(..., teststat="quadratic", 
distribution="asymptotic")}, as equivalents of \code{cctest(... ~ 1)} with 
chi-squared approximation and for handling data in blocks; implementations of 
canonical correlation analysis in other packages such as \pkg{CCP} (tests 
only), \pkg{MVar}, \pkg{candisc} (both including tests based on Wilks’ 
statistic), \pkg{yacca}, \pkg{CCA}, \pkg{acca}, \pkg{whitening}.
}

\examples{
## Artificial observations in 5-by-5 meter subplots in a forest for
## comparing cctest analyses with equivalent 'stats' methods:
dat <- within(data.frame(row.names=1:150), {
  u <- function() replicate(150, z<<-(z*69069+2^-32)\%\%1); z<-0
  plot <- factor(u() < .5, , c("a","b"))        # plot a or b
  x    <- as.integer(30*u() + c(1,82)[plot])    # x position on grid
  y    <- as.integer(30*u() + c(1,62)[plot])    # y position on grid
  ori  <- factor(u()\%/\%.25,,c("E","N","S","W")) # orientation of slope
  elev <- 40*u() + c(605,610)[plot]             # elevation (in meters)
  h    <- 115 - .15*elev + 2*log(1/u()-1)       # tree height (in meters)
  h5   <- h + log(1/u()-1)                      # tree height 5 years earlier
  h10  <- h5 + log(1/u()-1)                     # tree height 10 years earlier
  c15  <- as.integer(h10 + log(1/u()-1) > 20)   # 0-1 coded, 15 years earlier
  sapl <- as.integer(log(1/u())^.8*elev/40)     # number of saplings
  rm(u, z)
})

## t-tests:                                                   => p-value.beta
cctest(h~plot~1, dat)
  #t.test(h~plot, dat, var.equal=TRUE)
cctest(h-20~1~0, dat)
  #t.test(h~1, dat, mu=20)
cctest(h-h5~1~0, dat)
  #t.test(dat$h, dat$h5, paired=TRUE)

## Test for correlation:
cctest(h~elev~1, dat)
  #cor.test(~h+elev, dat)

## One-way analysis of variance:
cctest(h~ori~1, dat)
  #anova(lm(h~ori, dat))
% #oneway.test(h~ori, dat, var.equal=TRUE)  # some cases

## F-tests in linear models:
cctest(h~ori~1|elev, dat)
  #anova(lm(h~1+elev, dat), lm(h~ori+elev, dat))
cctest(h-h5~(h5-h10):(1|x|x^2)~0, dat)
  #summary(lm(h-h5~0+I(h5-h10)+I(h5-h10):poly(x,2,raw=TRUE), dat))

## Test in multivariate linear model based on Pillai's statistic:
cctest(h|h5|h10~x|y~1|elev, dat)
  #anova(lm(cbind(h,h5,h10)~elev, dat), lm(cbind(h,h5,h10)~elev+x+y, dat))

## Test based on Spearman's rank correlation coefficient:
cctest(rank(h)~rank(elev)~1, dat)
  #cor.test(~h+elev, dat, method="spearman", exact=FALSE)

## Kruskal-Wallis and Wilcoxon rank-sum tests:                => p-value.gamma
cctest(rank(h)~ori~1, dat)
  #kruskal.test(h~ori, dat)
cctest(rank(h)~plot~1, dat)
  #wilcox.test(h~plot, dat, exact=FALSE, correct=FALSE)

## Wilcoxon signed rank test:
cctest(rank(abs(h-h5))~sign(h-h5)~0, subset(dat, h-h5 != 0))
#dat|> within(d<-h-h5)|> subset(!!d)|> with(rank(abs(d))~sign(d)~0)|> cctest()
  #wilcox.test(h-h5 ~ 1, dat, exact=FALSE, correct=FALSE)

## Chi-squared test of independence:
cctest(ori~plot~1, dat, ~0)
cctest(ori~plot~1, as.data.frame(xtabs(~ori+plot,dat)), Freq~0)
  #summary(xtabs(~ori+plot, dat, drop.unused.levels=TRUE))
  #chisq.test(dat$ori, dat$plot, correct=FALSE)
% #prop.test(xtabs(~ori+plot,dat), correct=FALSE)  # some cases

## Score test in logistic regression (logit model, ...~1 only):
cctest(c15~x|y~1, dat, ~0)
  #anova(glm(c15~1, binomial, dat, epsilon=1e-20, maxit=200),
  #  glm(c15~1+x+y, binomial, dat), test="Rao")

## Score test in multinomial logit model (...~1 only):
cctest(ori~x|y~1, dat, ~0)
  #with(expand.grid(i=seq_len(nrow(dat)), j=unique(dat$ori)),
  #  anova(glm(ori==j ~ j+x+y, poisson, dat[i,], epsilon=1e-20, maxit=200),
  #    glm(ori==j ~ j*(x+y), poisson, dat[i,]), test="Rao"))

## Absolute values of (partial) correlation coefficients:     => estimate
cctest(h~elev~1, dat)
  #cor(dat$h, dat$elev)
cctest(h~elev~1|x|y, dat)
  #cov2cor(estVar(lm(cbind(h,elev)~1+x+y, dat)))
cctest(h~x|y|elev~1, dat)$estimate^2
  #summary(lm(h~1+x+y+elev, dat))$r.squared
%print(cctest(c(2,1)~diag(2)~0)$est, 18)  # inexactness warning (case u==t)

## Canonical correlations:
cctest(h|h5|h10~x|y~1, dat)
  #cancor(dat[c("x","y")],dat[c("h","h5","h10")])$cor

## Linear discriminant analysis:                              => various
with(cctest(h|h5|h10~ori~1, dat, ~ori), y / sqrt(1-estimate^2)[col(y)])[1:7,]
  #predict(MASS::lda(ori~h+h5+h10,dat))$x[1:7,]

## Correspondence analysis:
cctest(ori~plot~1, as.data.frame(xtabs(~ori+plot,dat)), Freq~0)[1:2]
  #MASS::corresp(~plot+ori, dat)

## Prediction in multivariate linear model:
with(cctest(h|h5|h10~1|x|y~0, dat, plot=="a"~0),
    x \%*\% diag(estimate,ncol(x),ncol(y)) \%*\% yinv)[1:7,]
  #predict(lm(cbind(h,h5,h10)~1+x+y, dat, subset=plot=="a"), dat)[1:7,]

## Other constructions:
cctest(ave(h,plot,FUN=rank)~ori~plot, aggregate(h~ori+plot,dat,mean))
  #friedman.test(h~ori|plot, aggregate(h~ori+plot,dat,mean))
cctest(ori:cut(elev,4)-cut(elev,4):ori~1~0, dat)
  #mcnemar.test(dat$ori, cut(dat$elev,4), correct=FALSE)
cctest(abs(rank(h)-mean(rank(h)))~plot~1, dat)
  #ansari.test(h~plot, dat, exact=FALSE)
cctest(ave((rank(h,,"first")-mean(rank(h)))^2,rank(h,,"min"))~plot~1, dat)
  #mood.test(h~plot, dat)
cctest(qnorm(rank(-abs(h-ave(h,ori,FUN=median)))/(length(h)+1)/2)~ori~1, dat)
  #fligner.test(h~ori, dat)
cctest(h~diag(plot=="a",length(plot))~plot, dat[1:7,])
  #var.test(h~plot, dat[1:7,], alternative="greater")
cctest(ori:sum(Freq)/Freq-1~1~0, as.data.frame(xtabs(~ori,dat)),
    (if(all(Freq)) Freq^2/sum(Freq)/c(.2,.3,.4,.1) else stop())~0)
  #chisq.test(xtabs(~ori,dat), p=c(.2,.3,.4,.1))
with(cctest({h|h5|h10;0;0;0} ~ {h|h5|h10;diag(3)} ~ {1|0*x;0;0;0},
           c(dat,`{`=rbind), ~ {h|h5|h10;diag(3)} | {1|0*x;0;0;0}),
    list(estimate/(s<-sqrt((1-estimate^2)*df.residual)), t(xinv*s)))
  #prcomp(~h+h5+h10, dat)

## Handling of additional arguments and edge cases:
cctest(1:150~x<6|y<6~1, c(dat,`:`=`:`,`|`=`|`), sapl/mean(sapl[sapl!=0])~1)
  #anova(lm(1:150~x<6|y<6, dat, weights=sapl))
cctest(h-h5~h10-h5~c15/c15, dat, (x<6)^NA~c15/c15)[1:2]
cctest(h~h10~offset(h5)+I(c15/c15), dat, , subset=x<6, NULL, stats=TRUE)[1:2]
  #scale(with(subset(dat,x<6), (cbind(h10,h)-h5)*c15/c15))
cctest(y~x~1, expand.grid(x=c(1,7)+.001, y=c(1,7)-.001), tol=0.6)[3:4]
  #lm(x+y~x+y, expand.grid(x=c(1,7)+.001, y=c(1,7)-.001), tol=0.6)
cctest(numeric()~NULL~NULL)
  #lm(numeric()~0)}
