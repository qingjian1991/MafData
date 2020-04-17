#' KiKs_cancer
#'
#' Analyses of selection using the dNdScv and dNdSloc models. Default parameters typically increase the performance of the method on cancer genomic studies. Reference files are currently only available for the GRCh37/hg19 version of the human genome.\cr
#' This is the function based on \code{\link[dndscv]{dndscv}} package. Here, we calculate the Ki/Ks(i=1:75) based on poisson regression model, similar with dNdSloc model in Amino acid levels. Except likelihood ratio test for selection, we also include a score test. see \strong{Details}.
#'
#' @author Original author is Inigo Martincorena (Wellcome Trust Sanger Institute)
#' @author Qingjian Chen Modified it and add Ki/Ks model. e-mail: \email{chenqingjian2010@@163.com}
#' @references Martincorena I, et al. (2017) Universal Patterns of Selection in Cancer and Somatic Tissues. Cell 171(5):1029-1041.
#' @references Greenman C, Wooster R, Futreal PA, Stratton MR, & Easton DF (2006) Statistical analysis of pathogenicity of somatic mutations in cancer. Genetics 173(4):2187-2198.
#'
#' @param mutations Environment input, Table of mutations (5 columns: sampleID, chr, pos, ref, alt). Only list independent events as mutations.
#' @param maf  file names, with columns: "chr","pos","ref","mut","sampleID", separated by "Tab", without header.
#' @param hg38to19 F; convert hg38 to 19; In default, mutations are based on hg19;
#' @param profile_mle profile estimated by mle of possion regression, the output \strong{mutrates} could as input of this params.
#' @param gene_list List of genes to restrict the analysis (use for targeted sequencing studies)
#' @param syntofew F (default). In case of too few mutations in Syn, thus we use the Ki1 to replace the SYN. This will only work when your data are filtered out of SYN mutations.
#' @param refdb Reference database (path to .rda file)
#' @param sm Substitution model (precomputed models are available in the data directory); sm = c("Type1_7_CpG","Type2_13_CpG","Type3_17","Type4_96","Type5_192"). see \strong{details} also.
#' @param kc List of a-priori known cancer genes (to be excluded from the indel background model)
#' @param cv Covariates (a matrix of covariates -columns- for each gene -rows-) [default: reference covariates] [cv=NULL runs dndscv without covariates]
#' @param max_muts_per_gene_per_sample If n<Inf, arbitrarily the first n mutations by chr position will be kept
#' @param max_coding_muts_per_sample Hypermutator samples often reduce power to detect selection
#' @param use_indel_sites Use unique indel sites instead of the total number of indels (it tends to be more robust)
#' @param min_indels Minimum number of indels required to run the indel recurrence module
#' @param maxcovs Maximum number of covariates that will be considered (additional columns in the matrix of covariates will be excluded)
#' @param constrain_wnon_wspl This constrains wnon==wspl (this typically leads to higher power to detect selection)
#' @param outp 1 default; Output: 1 = Global dN/dS values; 2 = Global dN/dS and dNdSloc, 3 =  Global dN/dS, dNdSloc and dNdScv, 4 = 3 + socre test . see \strong{details}.
#' @param KiKsmodel run Ki/Ks (default = T);
#'
#' @return 'dndscv' returns a list of objects:
#' @return Ki_mle: a list of variables. both likelihood ratio test and score test are applied.
#' @return globaldnds: globaldnds
#' @return sel_cv: dNdSloc
#' @return mutations: annotated mutations
#' @return p1: plot KI against UI.
#' @return UI.cor: correlation between 75-KI and UI.
#' @return mutrates:  mutation rate matrix, could as input of \strong{profile_mle}
#' @return mle_submodel: the estimated par from mle model of possion mutation rate regression
#' @return poissmodel: mle model of possion mutation rate regression
#' @return genemuts: mutations summarized by genes
#' @return par_KI: KI test by including w1-75 into one poisson regression model.
#' @return sel_loc: dNdScv
#' @return nbreg: negative binomial model for dNdScv
#'
#' @import rtracklayer
#' @import seqinr
#' @import Biostrings
#' @import MASS
#' @import GenomicRanges
#' @import IRanges
#' @import dndscv
#' @import ggplot2
#' @import ggpmisc
#' @import magrittr
#'
#' @details load data before running this code.
#' @details data("RefCDS_hg19_sum", package="MafData.Data")
#' @details data("AminoAcidDistance", package="MafData.Data")
#' @details \strong{score test} score test is more strict than likelihood test(LR).
#' @details \strong{sm}, Substitution model: 1) sm =1:5; sm = c("Type1_7_CpG","Type2_13_CpG","Type3_17","Type4_96","Type5_192") or 2) sm ="192r_3w", "12r_3w","2r_3w" or created by yourselves.
#' @export

KiKs_cancer = function(mutations = NULL , maf =NULL , hg38to19=F, profile_mle=NULL ,gene_list = NULL, syntofew = F,refdb = "hg19", sm = "192r_3w", kc = "cgc81", cv = "hg19", max_muts_per_gene_per_sample = 10, max_coding_muts_per_sample = 3000, use_indel_sites = T, min_indels = 5, maxcovs = 20, constrain_wnon_wspl = T, outp = 1, KiKsmodel= T) {


  # add confidence interval and p-value by score test.-------------------------------------
  add_interval_score_test= function(m,s, M, S, adjust_p_value=F, method="BH"){

    ##Calculate conf interval using Katz method
    p1<-(m/(M+1))
    p2<-(s/(S+1))
    globaldnds<-p1/p2
    N1 <- M
    N2 <- S

    SE = sqrt( (1-p1)/(N1*p1) + (1-p2)/(N2*p2) )

    finalLowCI = globaldnds * exp(-1.96*SE)
    finalHighCI = globaldnds * exp(1.96*SE)

    ## calculate p-value based on score test.
    t=(m+s)
    Ta=(M+S)
    U=m-(t*(M/Ta))
    V=t*(M*((Ta-M)/Ta^2))
    testscore=U^2/V
    pval_SSB=pchisq(testscore,df=1,lower.tail=FALSE)

    df.global<-as.data.frame(cbind(globaldnds,finalLowCI,finalHighCI, pval_SSB))
    colnames(df.global)<-c("globaldnds","low_CI","high_CI","pval_score")

    if(adjust_p_value){
      df.global$pval_score.adj<-p.adjust(df.global$pval, method = method, n=nrow(df.global))
    }

    return(df.global)
  }


  ## 1. Environment
  message("[1] Loading the environment...")

  # [Input] Reference database
  if (refdb == "hg19") {
    #data("RefCDS_hg19_sum", package="dndscv")
  } else {
    load(refdb)
  }
  data("AminoAcidDistance", package="MafData.Data")


  MutTable$N_EI = array(0, dim = c(192, 75))
  # [Input] Gene list (The user can input a gene list as a character vector)
  if (is.null(gene_list)) {
    gene_list = sapply(RefCDS, function(x) x$gene_name) # All genes [default]
  } else { # Using only genes in the input gene list
    allg = sapply(RefCDS,function(x) x$gene_name)
    nonex = gene_list[!(gene_list %in% allg)]
    if (length(nonex)>0) { message(sprintf("Warning
                                           The following input gene names are not in the RefCDS database: %s", paste(nonex,collapse=", ")))    }
    RefCDS = RefCDS[allg %in% gene_list] # Only input genes
    gr_genes = gr_genes[gr_genes$names %in% gene_list] # Only input genes
    message(sprintf("Length of RefCDS is %s",length(RefCDS)))
    ##### re-calculate L_EI
    Lall_EI = array(sapply(RefCDS, function(x) x$L_EI), dim=c(192,75,length(RefCDS)))
    L_EI = apply(Lall_EI, c(1,2), sum)
    MutTable$L_EI= L_EI
    }

  N_mut = matrix(data = 0, nrow = length(RefCDS), ncol = 75 )

  # [Input] Covariates (The user can input a custom set of covariates as a matrix)
  if (is.character(cv)) {
    data(list=sprintf("covariates_%s",cv), package="dndscv")
  } else {
    covs = cv
  }

  # [Input] Known cancer genes (The user can input a gene list as a character vector)
  if (kc[1] %in% c("cgc81")) {
    data(list=sprintf("cancergenes_%s",kc), package="dndscv")
  } else {
    known_cancergenes = kc
  }

  # [Input] Substitution model (The user can also input a custom substitution model as a matrix)

  makeSubstmodel = function(submod=submod, sm=1, colnm=NULL ){
    # trans Mutation_Type into sm.
    #ms = 1:5 or colnm for new column names.
    cols=c("Type1_7_CpG","Type2_13_CpG","Type3_17","Type4_96","Type5_192")
    if(is.null(colnm)){
      if(sm<=5){colnm = cols[sm]}
    }
    message(sprintf("sm is %s",colnm))
    mut = submod[,colnm]
    mut = paste("t*",mut, sep = "")
    chr =  mut[length(mut)]
    mut[mut %in% chr]= "t"; #change the last one into t.

    dt = data.frame(V2=mut,
                    V3=sprintf("%s*wmis",mut) ,
                    V4=sprintf("%s*wnon",mut) ,
                    V5=sprintf("%s*wspl",mut) )
    dt = as.matrix(dt)
  }

  if (length(sm)==1) {

    if(typeof(sm)=="character"){
      data(list=sprintf("submod_%s",sm), package="dndscv")
    }else if(typeof(sm)=="double"){
      data("submod", package = "MafData.Data")
      substmodel = makeSubstmodel(submod, sm=sm)
    }

  } else {
    substmodel = sm
  }


  ## 2. Mutation annotation
  message("[2] Annotating the mutations...")

  ####  reading  Mutations.

  if( !is.null(mutations) & is.null(maf) ) {
    #data("dataset_simbreast", package="dndscv");
    colnames(mutations) = c("sampleID","chr","pos","ref","mut")
  }else if(maf == "test"){
    data("dataset_simbreast", package="dndscv");
    colnames(mutations) = c("sampleID","chr","pos","ref","mut")

  }else{
    mutations=read.table(maf, sep="\t", header = F, stringsAsFactors = F)
    colnames(mutations) = c("chr","pos","ref","mut","sampleID")
    mutations=mutations[ ,c("sampleID","chr","pos","ref","mut")]
  }

  mutations[, c(1, 2, 4, 5)] = lapply(mutations[, c(1, 2, 4, 5)], as.character)




  #####
  if(hg38to19==T){
    ########################
    ###### Trans Mutations from hg38 to hg19.

    Tran38to19=function(mutations){
      gr_muts = GRanges(mutations$chr, IRanges(mutations$pos,mutations$pos), mcols=DataFrame(mutations[,c("sampleID","ref","mut")]) )
      genome(gr_muts)="hg19"
      data("chain", package="MafData.Data");
      #chain = import.chain("hg38ToHg19.over.chain")
      seqlevelsStyle(gr_muts) = "UCSC"
      gr_muts_ch = unlist( liftOver(gr_muts, chain38) )
      genome(gr_muts_ch)="hg38"

      mutations=data.frame(gr_muts_ch)[, c("mcols.sampleID","seqnames", "start", "mcols.ref", "mcols.mut")]
      colnames(mutations) = c("sampleID","chr","pos","ref","mut")
      mutations[,c(1,2,4,5)] = lapply(mutations[,c(1,2,4,5)], as.character) # Factors to character
      mutations$chr = sub("chr","",mutations$chr, fixed=T)
      return(mutations)
    }
    mutations=Tran38to19(mutations)
  }

  #### remove all indels.
  #nt = c("A","C","G","T")

  #nmutations = nrow(mutations)
  #mutations =dplyr::filter(mutations, ref %in% nt  & mut %in% nt )

  #message( sprintf("%s mutations (%.4f removed) remain among %s total mutations",  nrow(mutations) , (nmutations-nrow(mutations))/nmutations ,nmutations ))


  nt = c("A","C","G","T")
  trinucs = paste(rep(nt,each=16,times=1),rep(nt,each=4,times=4),rep(nt,each=1,times=16), sep="")
  trinucinds = setNames(1:64, trinucs)

  trinucsubs = NULL
  for (j in 1:length(trinucs)) {
    trinucsubs = c(trinucsubs, paste(trinucs[j], paste(substr(trinucs[j],1,1), setdiff(nt,substr(trinucs[j],2,2)), substr(trinucs[j],3,3), sep=""), sep=">"))
  }
  trinucsubsind = setNames(1:192, trinucsubs)

  ind = setNames(1:length(RefCDS), sapply(RefCDS,function(x) x$gene_name))
  gr_genes_ind = ind[gr_genes$names]

  # Warning about possible unannotated dinucleotide substitutions
  if (any(diff(mutations$pos)==1)) {
    warning("Mutations observed in contiguous sites within a sample. Please annotate or remove dinucleotide or complex substitutions for best results.")
  }

  # Warning about multiple instances of the same mutation in different sampleIDs
  if (nrow(unique(mutations[,2:5])) < nrow(mutations)) {
    warning("Same mutations observed in different sampleIDs. Please verify that these are independent events and remove duplicates otherwise.")
  }

  # Mapping mutations to genes
  gr_muts = GRanges(mutations$chr, IRanges(mutations$pos,mutations$pos))
  ol = as.matrix(findOverlaps(gr_muts, gr_genes, type="any", select="all"))
  mutations = mutations[ol[,1],] # Duplicating subs if they hit more than one gene
  mutations$geneind = gr_genes_ind[ol[,2]]
  mutations$gene = sapply(RefCDS,function(x) x$gene_name)[mutations$geneind]

  # Optional: Excluding samples exceeding the limit of mutations/sample [see Default parameters]
  nsampl = sort(table(mutations$sampleID))
  exclsamples = NULL
  if (any(nsampl>max_coding_muts_per_sample)) {
    message(sprintf('    Note: %0.0f samples excluded for exceeding the limit of mutations per sample',sum(nsampl>max_coding_muts_per_sample)))
    exclsamples = names(nsampl[nsampl>max_coding_muts_per_sample])
    mutations = mutations[!(mutations$sampleID %in% names(nsampl[nsampl>max_coding_muts_per_sample])),]
  }

  # Optional: Limiting the number of mutations per gene per sample (to minimise the impact of unannotated kataegis and other mutation clusters) [see Default parameters]
  mutrank = ave(mutations$pos, paste(mutations$sampleID,mutations$gene), FUN = function(x) rank(x))
  exclmuts = NULL
  if (any(mutrank>max_muts_per_gene_per_sample)) {
    message(sprintf('    Note: %0.0f mutations removed for exceeding the limit of mutations per gene per sample',sum(mutrank>max_muts_per_gene_per_sample)))
    exclmuts = mutations[mutrank>max_muts_per_gene_per_sample,]
    mutations = mutations[mutrank<=max_muts_per_gene_per_sample,]
  }

  # Additional annotation of substitutions

  snv = (mutations$ref %in% nt & mutations$mut %in% nt)
  indels = mutations[!snv,]
  mutations = mutations[snv,]
  mutations$ref_cod = mutations$ref
  mutations$mut_cod = mutations$mut
  compnt = setNames(rev(nt), nt)

  mutations$strand = sapply(RefCDS,function(x) x$strand)[mutations$geneind]
  isminus = (mutations$strand==-1)
  mutations$ref_cod[isminus] = compnt[mutations$ref[isminus]]
  mutations$mut_cod[isminus] = compnt[mutations$mut[isminus]]


  for (j in 1:length(RefCDS)) {
    RefCDS[[j]]$N = array(0, dim=c(192,4)) # Initialising the N matrices
  }

  # Subfunction: obtaining the codon positions of a coding mutation given the exon intervals

  chr2cds = function(pos,cds_int,strand) {
    if (strand==1) {
      return(which(pos==unlist(apply(cds_int, 1, function(x) x[1]:x[2]))))
    } else if (strand==-1) {
      return(which(pos==rev(unlist(apply(cds_int, 1, function(x) x[1]:x[2])))))
    }
  }

  # Annotating the functional impact of each substitution and populating the N matrices

  ref3_cod = mut3_cod = wrong_ref = aachange = ntchange = impact = old_cod= new_cod= array(NA, nrow(mutations))

  for (j in 1:nrow(mutations)) {

    geneind = mutations$geneind[j]
    pos = mutations$pos[j]
    if (any(pos == RefCDS[[geneind]]$intervals_splice)) { # Essential splice-site substitution

      impact[j] = "Essential_Splice"; impind = 4
      pos_ind = (pos==RefCDS[[geneind]]$intervals_splice)
      cdsnt = RefCDS[[geneind]]$seq_splice[pos_ind]
      ref3_cod[j] = sprintf("%s%s%s", RefCDS[[geneind]]$seq_splice1up[pos_ind], RefCDS[[geneind]]$seq_splice[pos_ind], RefCDS[[geneind]]$seq_splice1down[pos_ind])
      mut3_cod[j] = sprintf("%s%s%s", RefCDS[[geneind]]$seq_splice1up[pos_ind], mutations$mut_cod[j], RefCDS[[geneind]]$seq_splice1down[pos_ind])
      aachange[j] = ntchange[j] = "-"

    } else { # Coding substitution

      pos_ind = chr2cds(pos, RefCDS[[geneind]]$intervals_cds, RefCDS[[geneind]]$strand)
      cdsnt = RefCDS[[geneind]]$seq_cds[pos_ind]
      ref3_cod[j] = sprintf("%s%s%s", RefCDS[[geneind]]$seq_cds1up[pos_ind], RefCDS[[geneind]]$seq_cds[pos_ind], RefCDS[[geneind]]$seq_cds1down[pos_ind])
      mut3_cod[j] = sprintf("%s%s%s", RefCDS[[geneind]]$seq_cds1up[pos_ind], mutations$mut_cod[j], RefCDS[[geneind]]$seq_cds1down[pos_ind])
      codon_pos = c(ceiling(pos_ind/3)*3-2, ceiling(pos_ind/3)*3-1, ceiling(pos_ind/3)*3)
      old_codon = as.character(as.vector(RefCDS[[geneind]]$seq_cds[codon_pos]))
      pos_in_codon = pos_ind-(ceiling(pos_ind/3)-1)*3
      new_codon = old_codon;
      new_codon[pos_in_codon] = mutations$mut_cod[j]
      old_aa = seqinr::translate(old_codon)
      new_aa = seqinr::translate(new_codon)
      aachange[j] = sprintf('%s%0.0f%s',old_aa,ceiling(pos_ind/3),new_aa)
      ntchange[j] = sprintf('%s%0.0f%s',mutations$ref_cod[j],pos_ind,mutations$mut_cod[j])
      old_cod[j]= paste(old_codon, collapse = ""); new_cod[j]=paste(new_codon, collapse = "");
      # Annotating the impact of the mutation
      if (new_aa == old_aa){
        impact[j] = "Synonymous"; impind = 1
      } else if (new_aa == "*"){
        impact[j] = "Nonsense"; impind = 3
      } else if (old_aa != "*"){
        impact[j] = "Missense"; impind = 2

        ########################
        ####### In This Part, I add 75 aa.
        aa=sprintf("%s%s",new_aa, old_aa)
        aaid= ifelse(any(AminoAcidDistance$AA==aa), which(AminoAcidDistance$AA==aa),  which(AminoAcidDistance$AA==reverse(aa)) )

      } else if (old_aa=="*") {
        impact[j] = "Stop_loss"; impind = NA
      }
    }

    if (mutations$ref_cod[j] != as.character(cdsnt)) { # Incorrect base annotation in the input mutation file (the mutation will be excluded with a warning)
      wrong_ref[j] = 1
    } else if (!is.na(impind)) { # Correct base annotation in the input mutation file
      trisub = trinucsubsind[ paste(ref3_cod[j], mut3_cod[j], sep=">") ]
      RefCDS[[geneind]]$N[trisub,impind] = RefCDS[[geneind]]$N[trisub,impind] + 1 # Adding the mutation to the N matrices

      #add EIs
      if(impind == 2){
        MutTable$N_EI[trisub,aaid] = MutTable$N_EI[trisub,aaid] + 1 # Adding the mutation to the N matrices
        N_mut[geneind,aaid] = N_mut[geneind,aaid]+1
      }
    }

    if (round(j/1e4)==(j/1e4)) { message(sprintf('    %0.3g %%...', round(j/nrow(mutations),2)*100)) }
  }

  mutations$ref3_cod = ref3_cod
  mutations$mut3_cod = mut3_cod
  mutations$old_codon =old_cod
  mutations$new_codon = new_cod
  mutations$aachange = aachange
  mutations$ntchange = ntchange
  mutations$impact = impact
  mutations$pid = sapply(RefCDS,function(x) x$protein_id)[mutations$geneind]

  if (any(!is.na(wrong_ref))) {
    #stop(sprintf('%0.0f mutations have a wrong reference base, please correct and rerun.',sum(!is.na(wrong_ref)))) # This can be made into a mere warning and the rest of the code will work
    message(sprintf('Warning!!!
                    %0.0f mutations have a wrong reference base, please correct and rerun.',sum(!is.na(wrong_ref)))) # This can be made into a mere warning and the rest of the code will work

    wrong_refbase = mutations[!is.na(wrong_ref),]
    mutations = mutations[is.na(wrong_ref),]
  }

  if (any(nrow(indels))) { # If there are indels we concatenate the tables of subs and indels
    indels = cbind(indels, data.frame(ref_cod=".", mut_cod=".", strand=".", ref3_cod=".", mut3_cod=".", old_codon="." , new_codon =".", aachange=".", ntchange=".", impact="no-SNV", pid=sapply(RefCDS,function(x) x$protein_id)[indels$geneind]))
    annot = rbind(mutations, indels)
  } else {
    annot = mutations
  }
  annot = annot[order(annot$sampleID, annot$chr, annot$pos),]


## 3. Estimation of the global rate and selection parameters-------------------------------------------------------------------

  #### In this step, it estimate the global substitution rate.  r(i->j)
  message("[3] Estimating global rates...")

  # In this step, I suggest to reomve all known driver genes, thus get a real neutral mutation substitution profile. Because there are multi-hits in the same position of driver genes acorss different patients, it's apparent due to selection rather than mutation effect. So a real neutral profile should be removed all known driver genes.

  # Subfunction: fitting substitution model

  fit_substmodel = function(N, L, substmodel) {

    l = c(L); n = c(N); r = c(substmodel)
    n = n[l!=0]; r = r[l!=0]; l = l[l!=0]

    params = unique(base::strsplit(x=paste(r,collapse="*"), split="\\*")[[1]])
    indmat = as.data.frame(array(0, dim=c(length(r),length(params))))
    colnames(indmat) = params
    for (j in 1:length(r)) {
      indmat[j, base::strsplit(r[j], split="\\*")[[1]]] = 1
    }

    model = glm(formula = n ~ offset(log(l)) + . -1, data=indmat, family=poisson(link=log))
    mle = exp(coefficients(model)) # Maximum-likelihood estimates for the rate params
    ci = exp(confint.default(model)) # Wald confidence intervals
    par = data.frame(name=gsub("\`","",rownames(ci)), mle=mle[rownames(ci)], cilow=ci[,1], cihigh=ci[,2])
    return(list(par=par, model=model))
  }


  #get a neutral mutation rate or from input.
  if(is.null(profile_mle) ){

    # Removing drive genes and get a neutral mutation profile.
    allg = sapply(RefCDS,function(x) x$gene_name)
    RefCDS_neutral = RefCDS[ !(allg %in% known_cancergenes) ]
    Lall = array( sapply(RefCDS_neutral, function(x) x$L), dim=c(192,4,length(RefCDS_neutral)))
    Nall = array(sapply(RefCDS_neutral, function(x) x$N), dim=c(192,4,length(RefCDS_neutral)))
    L = apply(Lall, c(1,2), sum)
    N = apply(Nall, c(1,2), sum)

    poissout = fit_substmodel(N[,1:2], L[,1:2], substmodel[,1:2]) # Original substitution model #only SYN and nonSYN
    par = poissout$par
    poissmodel = poissout$model
    parmle =  setNames(par[,2], par[,1])
    mle_submodel = par
    rownames(mle_submodel) = NULL
    mutrates = sapply(substmodel[,1], function(x) prod(parmle[base::strsplit(x,split="\\*")[[1]]])) # Expected rate per available site
    numrates = length(mutrates)
  }else{
    message("inputting pre-calculated numrates")
    mutrates = profile_mle
    numrates = length(mutrates)
  }


  # Fitting all mutation rates and the 3 global selection parameters
  Lall = array(sapply(RefCDS, function(x) x$L), dim=c(192,4,length(RefCDS)))
  Nall = array(sapply(RefCDS, function(x) x$N), dim=c(192,4,length(RefCDS)))
  L = apply(Lall, c(1,2), sum)
  N = apply(Nall, c(1,2), sum)

  poissout = fit_substmodel(N, L, substmodel) # Original substitution model
  par = poissout$par
  poissmodel = poissout$model
  parmle =  setNames(par[,2], par[,1])
  mle_submodel = par
  rownames(mle_submodel) = NULL

  # Fitting models with 1 and 2 global selection parameters

  s1 = gsub("wmis","wall",gsub("wnon","wall",gsub("wspl","wall",substmodel)))
  par1 = fit_substmodel(N, L, s1)$par # Substitution model with 1 selection parameter
  s2 = gsub("wnon","wtru",gsub("wspl","wtru",substmodel))
  par2 = fit_substmodel(N, L, s2)$par # Substitution model with 1 selection parameter
  globaldnds = rbind(par, par1, par2)[c("wmis","wnon","wspl","wtru","wall"),]
  sel_loc = sel_cv = NULL

# -------------------------------------------------------------------

# 4. Running dNdSloc -------------------------------------------------------------------

  ## 4. dNdSloc: variable rate dN/dS model (gene mutation rate inferred from synonymous subs in the gene only)

  genemuts = data.frame(gene_name = sapply(RefCDS, function(x) x$gene_name), n_syn=NA, n_mis=NA, n_non=NA, n_spl=NA, exp_syn=NA, exp_mis=NA, exp_non=NA, exp_spl=NA)
  genemuts[,2:5] = t(sapply(RefCDS, function(x) colSums(x$N)))
  genemuts[,6:9] = t(sapply(RefCDS, function(x) colSums(x$L*mutrates)))

  if (outp > 1) {
    message("[4] Running dNdSloc...")

    selfun_loc = function(j) {
      y = as.numeric(genemuts[j,-1])
      x = RefCDS[[j]]

      # a. Neutral model: wmis==1, wnon==1, wspl==1
      mrfold = sum(y[1:4])/sum(y[5:8]) # Correction factor of "t" based on the obs/exp ratio of "neutral" mutations under the model
      ll0 = sum(dpois(x=x$N, lambda=x$L*mutrates*mrfold*t(array(c(1,1,1,1),dim=c(4,numrates))), log=T)) # loglik null model

      # b. Missense model: wmis==1, free wnon, free wspl
      mrfold = max(1e-10, sum(y[c(1,2)])/sum(y[c(5,6)])) # Correction factor of "t" based on the obs/exp ratio of "neutral" mutations under the model
      wfree = y[3:4]/y[7:8]/mrfold; wfree[y[3:4]==0] = 0
      llmis = sum(dpois(x=x$N, lambda=x$L*mutrates*mrfold*t(array(c(1,1,wfree),dim=c(4,numrates))), log=T)) # loglik free wmis

      # c. free wmis, wnon and wspl
      mrfold = max(1e-10, y[1]/y[5]) # Correction factor of "t"
      w = y[2:4]/y[6:8]/mrfold; w[y[2:4]==0] = 0 # MLE of dN/dS based on the local rate (using syn muts as neutral)
      llall = sum(dpois(x=x$N, lambda=x$L*mutrates*mrfold*t(array(c(1,w),dim=c(4,numrates))), log=T)) # loglik free wmis, wnon, wspl
      w[w>1e4] = 1e4

      p = 1-pchisq(2*(llall-c(llmis,ll0)),df=c(1,3))
      return(c(w,p))
    }

    sel_loc = as.data.frame(t(sapply(1:nrow(genemuts), selfun_loc)))
    colnames(sel_loc) = c("wmis_loc","wnon_loc","wspl_loc","pmis_loc","pall_loc")
    sel_loc$qmis_loc = p.adjust(sel_loc$pmis_loc, method="BH")
    sel_loc$qall_loc = p.adjust(sel_loc$pall_loc, method="BH")
    sel_loc = cbind(genemuts[,1:5],sel_loc)
    sel_loc = sel_loc[order(sel_loc$pall_loc,sel_loc$pmis_loc,-sel_loc$wmis_loc),]
  }

  # -------------------------------------------------------------------

#5. Running dNdSlcv-------------------------------------------------------------------

  ## 5. dNdScv: Negative binomial regression (with or without covariates) + local synonymous mutations

  nbreg = nbregind = NULL
  if (outp > 2) {

    message("[5] Running dNdScv...")

    # Covariates
    if (is.null(cv)) {
      nbrdf = genemuts[,c("n_syn","exp_syn")]
      model = MASS::glm.nb(n_syn ~ offset(log(exp_syn)) - 1 , data = nbrdf)
      message(sprintf("    Regression model for substitutions: no covariates were used (theta = %0.3g).", model$theta))
    } else {
      covs = covs[genemuts$gene_name,]
      if (ncol(covs) > maxcovs) {
        warning(sprintf("More than %s input covariates. Only the first %s will be considered.", maxcovs, maxcovs))
        covs = covs[,1:maxcovs]
      }
      nbrdf = cbind(genemuts[,c("n_syn","exp_syn")], covs)

      # Negative binomial regression
      model = suppressWarnings(MASS::glm.nb(n_syn ~ offset(log(exp_syn)) + . , data = nbrdf))
      if (!is.null(model$th.warn) | nrow(genemuts)<500) { # If there are warnings or if <500 genes, we run the regression without covariates
        model = MASS::glm.nb(n_syn ~ offset(log(exp_syn)) - 1 , data = nbrdf)
        message(sprintf("    Regression model for substitutions: no covariates were used (theta = %0.3g).", model$theta))
      } else {
        message(sprintf("    Regression model for substitutions: all covariates were used (theta = %0.3g).", model$theta))
      }
    }
    if (all(model$y==genemuts$n_syn)) {
      genemuts$exp_syn_cv = model$fitted.values
    }
    theta = model$theta
    nbreg = model

    # Subfunction: Analytical opt_t using only neutral subs
    mle_tcv = function(n_neutral, exp_rel_neutral, shape, scale) {
      tml = (n_neutral+shape-1)/(exp_rel_neutral+(1/scale))
      if (shape<=1) { # i.e. when theta<=1
        tml = max(shape*scale,tml) # i.e. tml is bounded to the mean of the gamma (i.e. y[9]) when theta<=1, since otherwise it takes meaningless values
      }
      return(tml)
    }

    # Subfunction: dNdScv per gene
    selfun_cv = function(j) {
      y = as.numeric(genemuts[j,-1])
      x = RefCDS[[j]]
      exp_rel = y[5:8]/y[5]
      # Gamma
      shape = theta
      scale = y[9]/theta

      # a. Neutral model
      indneut = 1:4 # vector of neutral mutation types under this model (1=synonymous, 2=missense, 3=nonsense, 4=essential_splice)
      opt_t = mle_tcv(n_neutral=sum(y[indneut]), exp_rel_neutral=sum(exp_rel[indneut]), shape=shape, scale=scale)
      mrfold = max(1e-10, opt_t/y[5]) # Correction factor of "t" based on the obs/exp ratio of "neutral" mutations under the model
      ll0 = sum(dpois(x=x$N, lambda=x$L*mutrates*mrfold*t(array(c(1,1,1,1),dim=c(4,numrates))), log=T)) + dgamma(opt_t, shape=shape, scale=scale, log=T) # loglik null model

      # b. Missense model: wmis==1, free wnon, free wspl
      indneut = 1:2
      opt_t = mle_tcv(n_neutral=sum(y[indneut]), exp_rel_neutral=sum(exp_rel[indneut]), shape=shape, scale=scale)
      mrfold = max(1e-10, opt_t/sum(y[5])) # Correction factor of "t" based on the obs/exp ratio of "neutral" mutations under the model
      wfree = y[3:4]/y[7:8]/mrfold; wfree[y[3:4]==0] = 0
      llmis = sum(dpois(x=x$N, lambda=x$L*mutrates*mrfold*t(array(c(1,1,wfree),dim=c(4,numrates))), log=T)) + dgamma(opt_t, shape=shape, scale=scale, log=T) # loglik free wmis

      # c. Truncating muts model: free wmis, wnon==wspl==1
      indneut = c(1,3,4)
      opt_t = mle_tcv(n_neutral=sum(y[indneut]), exp_rel_neutral=sum(exp_rel[indneut]), shape=shape, scale=scale)
      mrfold = max(1e-10, opt_t/sum(y[5])) # Correction factor of "t" based on the obs/exp ratio of "neutral" mutations under the model
      wfree = y[2]/y[6]/mrfold; wfree[y[2]==0] = 0
      lltrunc = sum(dpois(x=x$N, lambda=x$L*mutrates*mrfold*t(array(c(1,wfree,1,1),dim=c(4,numrates))), log=T)) + dgamma(opt_t, shape=shape, scale=scale, log=T) # loglik free wmis

      # d. Free selection model: free wmis, free wnon, free wspl
      indneut = 1
      opt_t = mle_tcv(n_neutral=sum(y[indneut]), exp_rel_neutral=sum(exp_rel[indneut]), shape=shape, scale=scale)
      mrfold = max(1e-10, opt_t/sum(y[5])) # Correction factor of "t" based on the obs/exp ratio of "neutral" mutations under the model
      wfree = y[2:4]/y[6:8]/mrfold; wfree[y[2:4]==0] = 0
      llall_unc = sum(dpois(x=x$N, lambda=x$L*mutrates*mrfold*t(array(c(1,wfree),dim=c(4,numrates))), log=T)) + dgamma(opt_t, shape=shape, scale=scale, log=T) # loglik free wmis

      if (constrain_wnon_wspl == 0) {

        p = 1-pchisq(2*(llall_unc-c(llmis,lltrunc,ll0)),df=c(1,2,3))
        return(c(wfree,p))

      } else { # d2. Free selection model: free wmis, free wnon==wspl

        wmisfree = y[2]/y[6]/mrfold; wmisfree[y[2]==0] = 0
        wtruncfree = sum(y[3:4])/sum(y[7:8])/mrfold; wtruncfree[sum(y[3:4])==0] = 0
        llall = sum(dpois(x=x$N, lambda=x$L*mutrates*mrfold*t(array(c(1,wmisfree,wtruncfree,wtruncfree),dim=c(4,numrates))), log=T)) + dgamma(opt_t, shape=shape, scale=scale, log=T) # loglik free wmis, free wnon==wspl
        p = 1-pchisq(2*c(llall_unc-llmis,llall-c(lltrunc,ll0)),df=c(1,1,2))
        return(c(opt_t, wmisfree,wtruncfree,wtruncfree,p))
      }
    }

    sel_cv = as.data.frame(t(sapply(1:nrow(genemuts), selfun_cv)))
    colnames(sel_cv) = c("opt_t","wmis_cv","wnon_cv","wspl_cv","pmis_cv","ptrunc_cv","pallsubs_cv")
    sel_cv$qmis_cv = p.adjust(sel_cv$pmis_cv, method="BH")
    sel_cv$qtrunc_cv = p.adjust(sel_cv$ptrunc_cv, method="BH")
    sel_cv$qallsubs_cv = p.adjust(sel_cv$pallsubs_cv, method="BH")
    sel_cv = cbind(genemuts[,1:10],sel_cv)
    #sel_cv = sel_cv[order(sel_cv$pallsubs_cv, sel_cv$pmis_cv, sel_cv$ptrunc_cv, -sel_cv$wmis_cv),] # Sorting genes in the output file

    ## Indel recurrence: based on a negative binomial regression (ideally fitted excluding major known driver genes)

    if (nrow(indels) >= min_indels) {

      geneindels = as.data.frame(array(0,dim=c(length(RefCDS),8)))
      colnames(geneindels) = c("gene_name","n_ind","n_induniq","n_indused","cds_length","excl","exp_unif","exp_indcv")
      geneindels$gene_name = sapply(RefCDS, function(x) x$gene_name)
      geneindels$n_ind = as.numeric(table(indels$gene)[geneindels[,1]]); geneindels[is.na(geneindels[,2]),2] = 0
      geneindels$n_induniq = as.numeric(table(unique(indels[,-1])$gene)[geneindels[,1]]); geneindels[is.na(geneindels[,3]),3] = 0

      if (use_indel_sites) {
        geneindels$n_indused = geneindels[,3]
      } else {
        geneindels$n_indused = geneindels[,2]
      }
      geneindels$cds_length = sapply(RefCDS, function(x) x$CDS_length)
      geneindels$excl = (geneindels[,1] %in% known_cancergenes)
      if (sum(geneindels[!geneindels$excl,"n_indused"]) == 0) { # If there are no indels for the background model we do not exclude any gene
        geneindels$excl = F
      }
      geneindels$exp_unif = sum(geneindels[!geneindels$excl,"n_indused"]) / sum(geneindels[!geneindels$excl,"cds_length"]) * geneindels$cds_length

      # Negative binomial regression for indels

      if (is.null(cv)) {
        nbrdf = geneindels[,c("n_indused","exp_unif")][!geneindels[,6],] # We exclude known drivers from the fit
        model = suppressWarnings(MASS::glm.nb(n_indused ~ offset(log(exp_unif)) + . , data = nbrdf))
        model = MASS::glm.nb(n_indused ~ offset(log(exp_unif)) - 1 , data = nbrdf)
        message(sprintf("    Regression model for indels: no covariates were used (theta = %0.3g)", model$theta))
        nbrdf_all = geneindels[,c("n_indused","exp_unif")]
      } else {
        nbrdf = cbind(geneindels[,c("n_indused","exp_unif")], covs)[!geneindels[,6],] # We exclude known drivers from the fit
        model = suppressWarnings(MASS::glm.nb(n_indused ~ offset(log(exp_unif)) + . , data = nbrdf))
        if (!is.null(model$th.warn) | nrow(genemuts)<500) { # If there are warnings or if <500 genes, we run the regression without covariates
          model = MASS::glm.nb(n_indused ~ offset(log(exp_unif)) - 1 , data = nbrdf)
          message(sprintf("    Regression model for indels: no covariates were used (theta = %0.3g)", model$theta))
        } else {
          message(sprintf("    Regression model for indels: all covariates were used (theta = %0.3g)", model$theta))
        }
        nbrdf_all = cbind(geneindels[,c("n_indused","exp_unif")], covs)
      }

      theta_indels = model$theta
      nbregind = model
      geneindels$exp_indcv = exp(predict(model,nbrdf_all))
      geneindels$wind = geneindels$n_indused / geneindels$exp_indcv

      # Statistical testing for indel recurrence per gene

      geneindels$pind = pnbinom(q=geneindels$n_indused-1, mu=geneindels$exp_indcv, size=theta_indels, lower.tail=F)
      geneindels$qind = p.adjust(geneindels$pind, method="BH")

      # Fisher combined p-values (substitutions and indels)

      sel_cv = merge(sel_cv, geneindels, by="gene_name")[,c("gene_name","n_syn","n_mis","n_non","n_spl","n_indused","wmis_cv","wnon_cv","wspl_cv","wind","pmis_cv","ptrunc_cv","pallsubs_cv","pind","qmis_cv","qtrunc_cv","qallsubs_cv")]
      colnames(sel_cv) = c("gene_name","n_syn","n_mis","n_non","n_spl","n_ind","wmis_cv","wnon_cv","wspl_cv","wind_cv","pmis_cv","ptrunc_cv","pallsubs_cv","pind_cv","qmis_cv","qtrunc_cv","qallsubs_cv")
      sel_cv$pglobal_cv = 1 - pchisq(-2 * (log(sel_cv$pallsubs_cv) + log(sel_cv$pind_cv)), df = 4)
      sel_cv$qglobal_cv = p.adjust(sel_cv$pglobal, method="BH")

      #sel_cv = sel_cv[order(sel_cv$pglobal_cv, sel_cv$pallsubs_cv, sel_cv$pmis_cv, sel_cv$ptrunc_cv, -sel_cv$wmis_cv),] # Sorting genes in the output file
    }
  }


  ############* add score test for each individual genes
  # need to do more. It is likeliy that score test is more strict than likelihood Ratio test(LR test)

  if(outp>3){
    message("[6] Running score test")

    #subfun for score test.
    average_mut_rate = sum(N[,1])/sum(L[,1])

    sel_cv = cbind(sel_cv,
                   add_interval_score_test(m =sel_cv$n_mis ,s =sel_cv$opt_t,
                                           M =sel_cv$exp_mis/average_mut_rate  ,S =sel_cv$exp_syn/average_mut_rate,
      adjust_p_value = T)

      )
    plot(sel_cv$pmis_cv, sel_cv$pval_score)
    abline(a=0,b=1,col="red")

  }


# -------------------------------------------------------------------


## 6. Estimation Ki/Ks from MLE methods.  -------------------------------------------------------------------

# In this step, the param **t** was assume constant in the genome.

  if(KiKsmodel ){
  message("[7] running Ki/Ks(i=1..75) by mle")

  N_EI = cbind(N[,1],  MutTable$N_EI );
  L_EI = cbind(L[,1],  MutTable$L_EI );
  kimuts = c( apply(N_EI, 2, sum), apply(L_EI*mutrates, 2, sum))

  #possion model, In this part, I add w1-75 to the general Poisson models. In this part, it use wald method to estimate params
  substmodel_EI =data.frame( V1 = substmodel[,1])
  for(i in 1:75 ){
    substmodel_EI = cbind(substmodel_EI, sprintf("%s*w%s",substmodel[,1], i) )
  }
  colnames(substmodel_EI) = sprintf("V%s",1:76)
  substmodel_EI= as.matrix( substmodel_EI)
  par_KI = fit_substmodel(N_EI, L_EI, substmodel_EI )$par %>% dplyr::filter(name %in% sprintf("w%s", 1:75) )


  selfun_Ki= function(Ki_rank=1, kimuts= kimuts, syntofew = F){
    if(syntofew){
      N_EI[,1]=N_EI[,2]
      L_EI[,1]=L_EI[,2]
      kimuts[1] = kimuts[2]
      kimuts[77] = kimuts[78]
    }
    #i= 1:75
    i = Ki_rank
    if(i<0 | i>75) stop("only 75 Ki's")
    average_mut_rate = sum(N_EI[,1])/sum(L_EI[,1])
    y = kimuts
    y1= y[1:76]  #obs_N
    y2 = y[77:152] #obs_L
    ## H0: W(i)==1; W(not i)!=1.
    mrfold = max(1e-10, sum(y1[c(1,i+1)])/sum(y2[c(1,i+1)])) # Correction factor of "t" based on the obs/exp ratio of "neutral" mutations under the model
    wfree = y1/y2/mrfold;
    wfree[c(1,i+1)]=1
    llmisi = sum(dpois(x=N_EI , lambda=L_EI*mutrates*mrfold*t(array( wfree ,dim=c(76,numrates))), log=T)) # loglik free wmis

    ## H1: free all.
    mrfold = max(1e-10, sum(y1[1])/sum(y2[1])) # Correction factor of "t" based on the obs/exp ratio of "neutral" mutations under the model
    wfree = y1/y2/mrfold;
    wfree[1]=1
    llall = sum(dpois(x=N_EI , lambda=L_EI*mutrates*mrfold*t(array( wfree ,dim=c(76,numrates))), log=T)) # loglik free wmis
    p = 1-pchisq(2*(llall- llmisi) ,df=1)

    return(c(y1[i+1],y2[i+1]/average_mut_rate,y1[1],mrfold*average_mut_rate, wfree[i+1],p))
  }

  Ki_mle = as.data.frame(t(sapply(1:75, selfun_Ki, kimuts= kimuts , syntofew=syntofew)))
  colnames(Ki_mle)=c("obs_N","exp_L","obs_S","Ks","KaKs","pval_mle")
  Ki_mle$qval_mle = p.adjust(Ki_mle$pval, method="BH")

  Ki_mle = cbind(Ki_mle,
                 add_interval_score_test(m=Ki_mle$obs_N, s=Ki_mle$obs_S,
                                         M=Ki_mle$exp_L, S=Ki_mle$obs_S/Ki_mle$Ks, adjust_p_value=T )
                  )

  if(syntofew){
   message("syntofew mode is active!!")
  }



# Summary: Plot and Correlation -------------------------------------------

  # This part is from "Ki_mle"

  ####################### Plot
  message("Summary: Plot and Correlation")

  AminoAcidDistance = AminoAcidDistance %>% mutate(Du = (max(UI)-UI)/(max(UI)-min(UI)) )
  dt.plot = data.frame(EI = Ki_mle$KaKs, Prop = AminoAcidDistance$Du)
  formula = 'y~x'

  UI.cor1 = cor.test(AminoAcidDistance$Du,Ki_mle$KaKs, method = "spearman")

  ypos = min( Ki_mle$KaKs) + 0.9* (max(Ki_mle$KaKs)- min(Ki_mle$KaKs))
  Ka.Ks =   globaldnds["wmis","mle"]

  formula= "y ~ x"
  p1 = ggplot( aes(x=Prop, y=EI ) , data=dt.plot) + theme_classic()  +
    geom_point(pch=21, fill="#3b518be5", size=2)+
    geom_smooth(method = "lm", formula = formula, se = F, col='blue') +
    geom_hline(aes(yintercept = Ka.Ks), linetype=2, col="red3", size=1 ) +
    annotate("text", x = 0.1 , y = ypos , label = paste("Rho is",
                                                        round( as.numeric(UI.cor1$estimate) , 3) )  ) +
    xlab(latex2exp::TeX("$\\Delta_U$"))+ylab("Observed Ki/Ks") +
    theme(text=element_text(size=14),
          axis.text=element_text(size=12),
          panel.spacing.x=unit(0, 'lines'),
          panel.spacing.y=unit(0, 'lines'),
          axis.line.x=element_line(color='darkgrey'),
          legend.position = "none"
    )

  ## estimate correlation between 75_KI and UI.
  UI.cor = cor.test(AminoAcidDistance$UI,Ki_mle$KaKs)
  UI.cor1 = cor.test(AminoAcidDistance$UI,Ki_mle$KaKs, method = "spearman")
  UI.cor=setNames(c(UI.cor$p.value, as.numeric(UI.cor$estimate ), UI.cor1$p.value, as.numeric(UI.cor1$estimate) ) , c("pval", "cor","pval1","rho"))


  }
# -------------------------------------------------------------------


if(outp>2){
        sel_cv = sel_cv[order(sel_cv$pallsubs_cv, sel_cv$pmis_cv, sel_cv$ptrunc_cv, -sel_cv$wmis_cv),] # Sorting genes in the output file
        if(nrow(indels) >= min_indels){
        sel_cv = sel_cv[order(sel_cv$pglobal_cv, sel_cv$pallsubs_cv, sel_cv$pmis_cv, sel_cv$ptrunc_cv, -sel_cv$wmis_cv),] # Sorting genes in the output file
        }

}




# -------------------------------------------


  #MutTable: Data for Count based Methods.

  dndscvout = list( globaldnds = globaldnds, mutations= mutations, mutrates=mutrates, mle_submodel= mle_submodel, poissmodel=poissmodel, genemuts =genemuts)

  if(outp >1 ){
    dndscvout$sel_loc=sel_loc
  }

  if(outp > 2){
    dndscvout$sel_cv=sel_cv
    dndscvout$nbreg =nbreg
  }

  if(KiKsmodel){
    dndscvout$Ki_mle=Ki_mle
    dndscvout$p1= p1
    dndscvout$UI.cor= UI.cor
    dndscvout$par_KI= par_KI
  }

  return(dndscvout)

}

