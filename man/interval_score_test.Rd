% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/interval_score_test.R
\name{interval_score_test}
\alias{interval_score_test}
\title{interval_score_test}
\usage{
interval_score_test(m, s, M, S, adjust_p_value = T, method = "BH")
}
\arguments{
\item{m, s, M, S}{m:missensen mutations; s: synonymous mutations; M: Sites for missense mutations ; S Sites for synonymous mutations.}

\item{adjust_p_value}{default(F); Whether calculates adjust p-values.}

\item{method}{method for adjust p-value. see \code{\link[stats]{p.adjust}}}
}
\value{
data.frame
}
\description{
Score test for Ka/Ks p-value estimation. Calculate confidence interval and p-value by score test.
}
\details{
p.adjust.methods include c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")
}
