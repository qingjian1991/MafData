#' interval_score_test
#'
#' Score test for Ka/Ks p-value estimation. Calculate confidence interval and p-value by score test.
#'
#' @param m,s,M,S m:missensen mutations; s: synonymous mutations; M: Sites for missense mutations ; S Sites for synonymous mutations.
#' @param adjust_p_value default(F); Whether calculates adjust p-values.
#' @param method method for adjust p-value. see \code{\link[stats]{p.adjust}}
#' @return data.frame
#' @details p.adjust.methods include c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")
#' @export


# calculate confidence interval and p-value by score test.-------------------------------------
interval_score_test= function(m, s, M, S, adjust_p_value=T, method="BH"){

  #chang Ka/Ks calculations.
  #pseudocount to add to keep neutrality =1
  pseudna<-M/S
  pn <- (m+pseudna)/M
  ps <- (s+0.5)/S
  estdnds <- pn/ps

  ##Calculate conf interval using Katz method
  p1<-(m+0.5)/(M+1)
  p2<-(s+0.5)/(S+1)
  globaldnds<-p1/p2
  N1 <- M
  N2 <- S

  SE = sqrt( (1-p1)/(N1*p1) + (1-p2)/(N2*p2) )

  finalLowCI = globaldnds * exp(-1.96*SE)
  finalHighCI = globaldnds * exp(1.96*SE)

  ## calculate p-value based on score test.
  #Define values of expected and observed and transform into a chi-square distribution
  t=(m+s)
  Ta=(M+S)
  U=m-(t*(M/Ta))
  V=t*(M*((Ta-M)/Ta^2))
  testscore=U^2/V
  pval_SSB=pchisq(testscore,df=1,lower.tail=FALSE)

  df.global<-as.data.frame(cbind(estdnds,globaldnds,finalLowCI,finalHighCI, pval_SSB))
  colnames(df.global)<-c("estdnds","globaldnds","low_CI","high_CI","pval_score")

  if(adjust_p_value){
    df.global$pval_score.adj<-p.adjust(df.global$pval, method = method, n=nrow(df.global))
  }

  return(df.global)
}


