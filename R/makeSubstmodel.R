#' makeSubstmodel
#'
#' create the new sm for dndscv sm models.
#'
#' @param submod load from the data: data("submod", package = "MafData")
#' @param sm \strong{sm}, Substitution model: 1) sm =1:5; sm = c("Type1_7_CpG","Type2_13_CpG","Type3_17","Type4_96","Type5_192")

#' @return substmodel for sm

#' @details This is for make new substitution models for dndscv.
#' @details you can create a new sm based on submod.

#' @export

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


