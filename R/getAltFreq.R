#' getAltFreq
#'
#' calculating the alternative allele frequency of mutations.
#'
#' @param MAF MAF file, with relevant fields consisten with the GDC MAF specification (\url{https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/}).
#' @param t_depth The colnames of t_depth.
#' @param t_alt_count The colnames of t_alt_count.
#' @return MAF after filtering.
#' @export

getAltFreq = function( MAF, t_depth = "t_depth", t_alt_count = "t_alt_count" ){

  if( !(t_depth %in%  colnames(MAF) ) ){
    stop( sprintf("There is no %s in MAF", t_depth))
  }

  if( !(t_alt_count %in%  colnames(MAF) ) ){
    stop( sprintf("There is no %s in MAF", t_alt_count))
  }

  MAF$t_alt_freq =mapply( mean, mapply(function(x1, x2){x1/x2}, lapply( strsplit(MAF[ ,t_alt_count] , split = "|", fixed=T), as.numeric ), lapply( strsplit(MAF[ ,t_depth] , split = "|", fixed=T), as.numeric ) ))

  return(MAF)

}















