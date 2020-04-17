#' getObsEI
#'
#'summary a lot of KiKs_cancer resluts into a list.
#'
#' @param list lists of KiKs_cancers
#' @param CancerName The corresponding cancer names
#'
#' @return lists including obsEI and obsGlobal
#' @export

getObsEI = function(list= list(), CancerName= NULL ){

  if(typeof(list)!="list"){ stop("data type is not list") }

  if(is.null(CancerName)){  stop("Cancer Names!!") }

  len= length(list)
  if(len == 0){stop("length of list is 0")  }

  if(len != length(CancerName)){ stop("the length of list and CancerName are not equal")}

  data =c()
  global = c()
  for(i in 1:len){

    data = rbind( data, data.frame(Cancer = CancerName[i],
                                   Distance = AminoAcidDistance$UI,
                                   value= list[[i]]$Ki_mle$KaKs,
                                   Ka.Ks =  list[[i]]$globaldnds["wmis","mle"],
                                   obsN= sum(list[[i]]$Ki_mle$obs_N)
    ))

    global = rbind(global, list[[i]]$globaldnds[c("wmis","wnon"),] %>%
                     mutate(Cancer = CancerName[i]) %>%
                     mutate(Cor = list[[i]]$UI.cor[4]) %>%
                     mutate(pval = list[[i]]$UI.cor[3] )
    )
  }
  data$label = sprintf( "%s\nn= %s; R= %.3f",data$Cancer, prettyNum( data$obsN ,big.mark = ','), data$Ka.Ks)
  data$label = factor(data$label, levels = unique(data$label))
  return(list(obsEI= data, obsGlobal= global)  )

}
