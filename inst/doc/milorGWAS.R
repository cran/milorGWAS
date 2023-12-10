## ----echo = FALSE-----------------------------------------------------------------------------------------------------
oldoptions <- options()
oldpar <- par()
options(width = 120)

## ----example11, echo=TRUE, message=FALSE------------------------------------------------------------------------------
library(milorGWAS)
x <- as.bed.matrix(TTN.gen, TTN.fam, TTN.bim)
x
options(gaston.verbose = FALSE)

## ----example12, echo=TRUE, message=FALSE------------------------------------------------------------------------------
set.seed(1)
## some covariates : an intercept, and a uniformly distributed covariate.
X <- cbind(1, runif(nrow(x)))
## A random GRM
ran <- random.pm(nrow(x))
## random effects (with variance tau = 1)
omega <- lmm.simu(1, 0, eigenK=ran$eigen)$omega
## linear term of the model
L <- X %*% c(0.1,-0.2) + omega
## vector of probabilities p = expit(L)
p <- 1/(1+exp( -L ))
## vector of binary phenotypes
y <- rbinom(length(p), 1, p)

## ----assoGaston, cache=TRUE-------------------------------------------------------------------------------------------
a.pql <- association.test(x, y, X, K = ran$K, method = "lmm", response = "bin", test = "wald")
a.gmm <- association.test(x, y, X, K = ran$K, method = "lmm", response = "bin", test = "score")

## ---------------------------------------------------------------------------------------------------------------------
head(a.pql)
head(a.gmm)

## ----assoMilor, cache=TRUE--------------------------------------------------------------------------------------------
a.ofst <- association.test.logistic(x, y, X, K = ran$K, algorithm = "offset")
a.amle <- association.test.logistic(x, y, X, K = ran$K, algorithm = "amle")
head(a.ofst)
head(a.amle)

## ----plotbeta, fig.height = 4, fig.width = 7, dev = 'png', dpi = 300--------------------------------------------------
par(mfrow=c(1,2), cex = 0.9)
plot(a.amle$beta, a.pql$beta, xlab = "AMLE", ylab = "PQL"); abline(0,1,col=4)
plot(a.ofst$beta, a.pql$beta, xlab = "Offset", ylab = "PQL"); abline(0,1,col=4)

## ----plotp, fig.height = 4, fig.width = 7, dev = 'png', dpi = 300-----------------------------------------------------
par(mfrow=c(1,2), cex = 0.9)
plot(a.amle$p, a.pql$p, log = "xy", xlab = "AMLE/score", ylab = "PQL"); abline(0,1,col=4)
plot(a.ofst$p, a.pql$p, log = "xy", xlab = "Offset", ylab = "PQL"); abline(0,1,col=4)

## ----eval=FALSE-------------------------------------------------------------------------------------------------------
#  install.packages("GridData", repos="https://genostats.github.io/R/")

## ----eval=FALSE-------------------------------------------------------------------------------------------------------
#  filepath <-system.file("extdata", "GridData.bed", package="GridData")
#  x <- read.bed.matrix(filepath)
#  x <- set.stats(x)

## ----eval=FALSE-------------------------------------------------------------------------------------------------------
#  x

## ----eval=FALSE-------------------------------------------------------------------------------------------------------
#  table(x@ped$pop)

## ----eval=FALSE-------------------------------------------------------------------------------------------------------
#  table(x@ped$pheno)

## ----eval=FALSE-------------------------------------------------------------------------------------------------------
#  table(x@ped$pop, x@ped$pheno)

## ----eval=FALSE-------------------------------------------------------------------------------------------------------
#  K <- GRM(x)
#  eigenK <- eigen(K)

## ----eval=FALSE-------------------------------------------------------------------------------------------------------
#  MLM <- association.test(x, method = "lmm", response = "quantitative",
#                          test = "wald", eigenK = eigenK, p = 10)
#  MLR <- association.test.logistic(x,  K = K, eigenK = eigenK, p = 10, algorithm = "amle")

## ----eval=FALSE-------------------------------------------------------------------------------------------------------
#  # SNPs categories from true strata information
#  cat.str <- SNP.category(x,x@ped$pop)
#  # SNPs categories from the first PCs coordinates
#  cat.PC1 <- SNP.category(x,-eigenK$vectors[,1])
#  cat.PC2 <- SNP.category(x,-eigenK$vectors[,2])

## ----eval=FALSE-------------------------------------------------------------------------------------------------------
#  par(mfrow=c(1,3), cex = 0.7)
#  qqplot.pvalues(MLM, cat.str, col.abline="blue", main="MLM - from strata", pch = 16)
#  qqplot.pvalues(MLM, cat.PC1, col.abline="blue", main="MLM - from PC1", pch = 16)
#  qqplot.pvalues(MLM, cat.PC2, col.abline="blue", main="MLM - from PC2", pch = 16)

## ----echo=FALSE-------------------------------------------------------------------------------------------------------
knitr::include_graphics("qqplotsMLM-1.png", dpi = 300)

## ----eval=FALSE-------------------------------------------------------------------------------------------------------
#  par(mfrow=c(1,3), cex = 0.7)
#  qqplot.pvalues(MLR, cat.str, col.abline="blue", main="MLR - from strata", pch = 16)
#  qqplot.pvalues(MLR, cat.PC1, col.abline="blue", main="MLR - from PC1", pch = 16)
#  qqplot.pvalues(MLR, cat.PC1, col.abline="blue", main="MLR - from PC2", pch = 16)

## ----echo=FALSE-------------------------------------------------------------------------------------------------------
knitr::include_graphics("qqplotsMLR-1.png", dpi = 300)

## ----echo = FALSE-------------------------------------------------------------
par(oldpar)
options(oldoptions)

