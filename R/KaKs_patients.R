#' KaKs_patients
#'
#' calculate the Ka/Ks per cancer case. The default is to use all genes across all patients to calculate the Ka/Ks per cancer case.
#'
#' @param dndscvOUT output from the dndscv, see \code{\link[dndscv]{dndscv}}
#' @param patientsID the patients used for the Ka/Ks calculation.
#' @param geneID the genes used for the Ka/Ks calculation.
#' @param sm Please set the sm same with pre-calculation.
#' @import tidyverse
#' @return KaKs_patients returns a list of objects:
#' @return - PatientsKaKs: patients Ka/Ks.
#' @return - globaldnds: globaldnds information.
#' @details The first step is calculate the calculate \code{\link[dndscv]{dndscv}} from the whole dataset, and then extract the information from the dndscv output.
#' @details The Ka/Ks per patients is equal to the number of A/S changes per the site of #A/#S. The A/S changes is from each patient and the #A/#S site is designated for the genome-wide #A/#S site.
#' @details The Ka/Ks can be calculated from a specific gene list rather a genome-wide Ka/Ks, thus we can set the geneID for the specific gene list.
#' @details In the other situation, the Ka/Ks is based on a certain patients rather than whole patients. This will affect the genome-wide #A/#S site. Set the patientsID for a certain group of patients.
#' @export

KaKs_patients=function(dndscvOUT, patientsID = NULL , geneID = NULL, sm = "192r_3w" ){


  #calculate the mutrates
  # [Input] Substitution model (The user can also input a custom substitution model as a matrix)
  #mutrates = dndscvOUT$mutrates
  if (length(sm)==1) {
    data(list=sprintf("submod_%s",sm), package="dndscv")
  } else {
    substmodel = sm
  }
  mle_submodel=dndscvOUT$mle_submodel
  parmle =  setNames(mle_submodel[,2], mle_submodel[,1])
  # Expected rate per available site
  mutrates = sapply(substmodel[,1], function(x) prod(parmle[base::strsplit(x,split="\\*")[[1]]]))


  #estimate the A/S per genes.
  annotmuts = dndscvOUT$annotmuts
  Estimate_AS_ratio = function(x){
    apply(x$L*mutrates, 2, sum)
  }
  Gene_AS_ratio = as.data.frame(t(sapply(RefCDS, Estimate_AS_ratio)))
  colnames(Gene_AS_ratio) = c("Lsyn","Lmis","Lnon","Lspl")
  Gene_AS_ratio$gene = sapply(RefCDS, function(x) x$gene_name )

  #table(dndscvOUT$annotmuts$impact)
  #Essential_Splice      Missense    Nonsense       Synonymous

  if(!is.null(geneID)){
    message("Use selected Genes to calculate the Ka/Ks per patient")
    annotmuts = annotmuts %>% filter(gene %in% geneID)
    Gene_AS_ratio = Gene_AS_ratio %>% filter(gene %in% geneID)
    allg = sapply(RefCDS,function(x) x$gene_name)
    RefCDS =RefCDS[allg %in% geneID]
  }

  if(!is.null(patientsID)){
    message("Use selected Patientes to calculate the global dnds")
    annotmuts = annotmuts %>% filter(sampleID  %in% sampleID)
  }

  #mu = n/L =Ks of synonymous.

  # L is restricted in the selected region.
  L = array(sapply(RefCDS, function(x) x$L), dim=c(192,4,length(RefCDS) )) %>%
    apply(c(1,2), sum)
  #average_mut_rate is Ks.
  average_mut_rate = annotmuts %>% filter(impact == "Synonymous") %>%
    nrow()/sum(L[,1])

  globaldnds = as.data.frame.array( t( apply(Gene_AS_ratio[,c(1:4)], 2, sum)))

  globaldnds = annotmuts %>%
    mutate(impact = factor(impact, levels = c("Synonymous","Missense","Nonsense","Essential_Splice")) ) %>%
    group_by(impact) %>%
    summarise(MutNum =  n()) %>%
    spread(key = impact, value = MutNum, fill = 0) %>%
    cbind(., globaldnds) %>%
    mutate(mu = average_mut_rate,
           wmis = ( Missense/Synonymous)/(Lmis/Lsyn),
           wnon = ( Nonsense/Synonymous)/(Lnon/Lsyn),
           wspl =( Essential_Splice/Synonymous)/(Lspl/Lsyn),
           wtru =( (Nonsense+Essential_Splice)/Synonymous)/((Lnon+Lspl)/Lsyn),
           wall =( (Nonsense+Essential_Splice+Missense)/Synonymous)/((Lnon+Lspl+Lmis)/Lsyn),
             )

  PatientsKaKs = annotmuts %>%
    mutate(impact = factor(impact, levels = c("Synonymous","Missense","Nonsense","Essential_Splice")) ) %>%
    group_by(sampleID, impact ) %>%
    summarise(MutNum =  n()) %>%
    spread(key = impact, value = MutNum, fill = 0) %>%
    mutate(
      wmis = ( Missense/(Synonymous+0.5))/(globaldnds$Lmis/globaldnds$Lsyn),
      wnon = ( Nonsense/(Synonymous+0.5))/(globaldnds$Lnon/globaldnds$Lsyn),
      wspl =( Essential_Splice/(Synonymous+0.5))/(globaldnds$Lspl/globaldnds$Lsyn),
      wtru =( (Nonsense+Essential_Splice)/(Synonymous+0.5))/((globaldnds$Lnon+globaldnds$Lspl)/globaldnds$Lsyn),
      wall =( (Nonsense+Essential_Splice+Missense)/(Synonymous+0.5))/((globaldnds$Lnon+globaldnds$Lspl+globaldnds$Lmis)/globaldnds$Lsyn)
    )

  #add p-values
  PatientsKaKs = interval_score_test(m = PatientsKaKs$Missense,
                      s = PatientsKaKs$Synonymous,
                      M = globaldnds$Lmis/globaldnds$mu,
                      S = globaldnds$Lsyn/globaldnds$mu
                      ) %>%
  cbind(data.frame(PatientsKaKs), .)

 return(list(
   PatientsKaKs = PatientsKaKs,
   globaldnds = globaldnds
 ))

}

