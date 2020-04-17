#' output formate
#'
#' Input is MAF, output is \tabular{lcr}{\cr}
#' 1) vcf; files, vcf  \tabular{lcr}{\cr}
#' 2) mut5.file; Files, simple 5 column, "chr","pos","ref","mut","sampleID", separated by \emph{"Tab"}, without header. \tabular{lcr}{\cr}
#' 3) mut5; system variables, simple 5 column. with header is "sampleID chr pos ref mut \tabular{lcr}{\cr}
#'
#' @param MAF MAF files
#' @param columns 5 column names for output(chr, pos, ref, mut, sampleID).
#' @param vcf  1) vcf file names
#' @param mut5.file 2) file names
#' @param mut5 3) T or F
#' @param cancer A simple variable, output cancer name in vcf output.
#' @param columns_Cancer  Column names in MAF for cancer type. The \emph{columns_Cancer} has \strong{higher priority} than \emph{cancer}.
#' @return mut5 and mut5 file are for dndscv. and vcf is for SSB_selection.
#' @export
#' @details \strong{cancer} and \strong{columns_Cancer}: Set Cancer Names in VCF. The Names could be from MAF (\strong{columns_Cancer}) or params  \strong{cancer}


out_formate= function(MAF,
                      columns= c("Chromosome","Start_Position","Reference_Allele","Tumor_Seq_Allele2","Tumor_Sample_Barcode"),
                      mut5=T, vcf = NULL, mut5.file= NULL , cancer="cancer", columns_Cancer="Cancer"){

  if(!is.null( mut5.file)) {
    mutations = MAF[,columns]
    message( sprintf( "output mut5 file: %s", mut5.file ) )
    write.table(mutations, file= mut5.file , row.names = F, col.names = F, quote = F, sep = "\t")
  }

  if(!is.null(vcf)){
    message( sprintf( "output vcf file: %s", vcf ) )

  #set Cancer Names in VCF, The Names could be from MAF (provided by columns_Cancer) or params cancer
  CancerName= cancer
  if( is.null(columns_Cancer) & cancer == "cancer"){
    message("Please changes cancer names in the vcf files if necessary!")
  }

  if(! is.null(columns_Cancer) ){
    CancerName= MAF[,columns_Cancer]
  }


    mut.vcf = MAF[, columns[1:2] ]
    mut.vcf$ID = "."
    mut.vcf = cbind( mut.vcf, MAF[,columns[3:4] ] )
    mut.vcf$QUAL = "."
    mut.vcf$FILTER = "."
    mut.vcf$INFO = sprintf( "ID=%s;CA=%s", MAF[,columns[5]], CancerName)

    colnames(mut.vcf)=c("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO")
    write.table(mut.vcf, file= vcf , row.names = F, col.names = T, quote = F, sep = "\t")

  }

  if(mut5){
    mutations = MAF[,columns]
    colnames(mutations)=  c("chr","pos","ref","mut","sampleID")
    mutations=mutations[ ,c("sampleID","chr","pos","ref","mut")]
    return(mutations)
  }

}

