% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/qqplot_pvalues.r
\name{qqplot.pvalues}
\alias{qqplot.pvalues}
\title{Stratified QQ-plot of p-values}
\usage{
qqplot.pvalues(
  p,
  snp.cat,
  col.cat,
  col.abline = "red",
  CB = TRUE,
  col.CB = "gray80",
  CB.level = 0.95,
  thinning = TRUE,
  ...
)
}
\arguments{
\item{p}{vector of p-values, or a data.frame with a column named \code{p}}

\item{snp.cat}{(optional) A factor giving the SNP categories.}

\item{col.cat}{(optional) A vector of colors used to plot the SNP categories.}

\item{col.abline}{Color of the line of slope 1. Set to \code{NA} to suppress.}

\item{CB}{\code{Logical}. If \code{TRUE}, a confidence band is included in the plot.}

\item{col.CB}{The color of the confidence band.}

\item{CB.level}{The level of the confidence band.}

\item{thinning}{\code{Logical}. If \code{TRUE}, not all points are displayed.}

\item{...}{Graphical parameters to be passed to \code{plot} and \code{points}}
}
\value{
Returns a 'NULL'
}
\description{
Draws a QQ plot of p-values
}
\details{
This function draws a QQ plot of \eqn{p}-values, stratified by categories.
If the parameter \code{snp.cat} is missing, the function falls back on \code{gaston::qqplot.pvalues}.
}
\examples{
# a random vector of categories
ca <- sample(c("A","B","C"), 1e6, TRUE, c(0.05, 0.9, 0.05))
# a vector of p-values, with different distribution depending on the strata
p <- runif(1e6)**ifelse(ca == "A", .8, ifelse(ca == "B", 1, 1.2))
qqplot.pvalues(p, ca)

}
\seealso{
\code{\link{SNP.category}}, \code{\link[gaston]{qqplot.pvalues}} (in gaston)
}
