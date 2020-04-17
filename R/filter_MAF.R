#' filter MAF
#'
#' Input is MAF, output is MAF.
#'
#' @param MAF MAF file, with relevant fields consisten with the GDC MAF specification (\url{https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/}).
#' @param filterMC3 whether filtering mutation not in MC3, default is FALSE. see details.
#' @param MC3_name The colnames of MC3. If the MC3 is TRUE, this mutation is retained.
#' @param addAltFreq whether add the alt frequency of mutations, default is FALSE.
#' @param t_depth The colnames of t_depth.
#' @param t_alt_count The colnames of t_alt_count.
#' @param maxNum patients with mutations more than this number will be removed.
#' @param minNum patients with mutations less than this number will be removed.
#' @param patientsID removing patients without in patients ID.
#' @param convert_coordinate  whether convert hg38 to 19(or hg19 to hg38), default is FALSE.
#' @param chain chain file "hg38ToHg19.over.chain"
#' @return MAF after filtering.
#' @details we convert data from hg38 to hg19 , keep only SNP data , removing hyper-mutated patients.
#' @details the functions are including: (1)hg_converter (2)unique_tumor_addition_function (3)tumor_allele_adder (4)DNP_TNP_remover (5)removing_patients (6)getAltFreq
#' @details \strong{MC3:The Multi-Center Mutation Calling in Multiple Cancers}
#' @details The MC3 data set provides consistent variant calling and filtering across the 10K patients in The Cancer Genome Atlas (TCGA). see \url{https://gdc.cancer.gov/about-data/publications/mc3-2017}
#' @export
#' @examples
#' maf <- read.delim(
#'  file = system.file("extdata", "TCGA.KIRC.merge.combine.uniq.amaf", package = "MafData") ,
#'  header = T,skip = 1,stringsAsFactors = F)
#' chain = system.file("extdata", "hg38ToHg19.over.chain", package = "MafData")
#' maf_left = filter_MAF(MAF = maf, hg38to19 = T, chain = chain )

filter_MAF = function(MAF, filterMC3 = F, MC3_name = "MC3_Overlap" ,addAltFreq = F, t_depth = "t_depth", t_alt_count = "t_alt_count", maxNum = 3000, minNum = 0,  patientsID= NULL,  convert_coordinate = F, chain="hg38ToHg19.over.chain") {

if(convert_coordinate){
message("1) hg_converter; ", chain)
MAF <- hg_converter(chain = chain, maf_to_convert = MAF)
}

message("2) unique_tumor_addition_function ;summary for tumor data")
MAF <- unique_tumor_addition_function(MAF.file = MAF)

message("3) tumor_allele_adder; check out tumor column names.")
MAF <- tumor_allele_adder(MAF = MAF)
print( sprintf("3) nrow of MAF: %s ", nrow(MAF)))

message( "4) DNP_TNP_remover; Removing DNV and TNV" )
MAF <- DNP_TNP_remover(MAF = MAF)
print( sprintf("4) nrow of MAF: %s ", nrow(MAF)))

message( "5)removing_patients;  removing patients.
# 1)removing patients.  2) removing patients more than N mutations.")

MAF <- removing_patients(MAF= MAF, patientsID = patientsID, maxNum = maxNum, minNum = minNum )

if(filterMC3){
  message("6) filterMC3")

  if(is.character(MAF[,MC3_name])){
    MAF[,MC3_name] = ifelse( toupper(MAF[,MC3_name]) == "TRUE", TRUE, FALSE )
  }

  if(is.logical(MAF[,MC3_name])){
    MAF <- MAF[ MAF[,MC3_name] == TRUE , ]
  }else{
    message( sprintf("The column %s is not TRUE or FALSE", MC3_name ))
  }

}

if(addAltFreq){
  message("7) getAltFreq")
  MAF <- getAltFreq(MAF = MAF, t_depth = t_depth, t_alt_count = t_alt_count)
}

return(MAF)

}




