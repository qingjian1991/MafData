#' removing patients
#'
#' removing patients belongs to hyper-mutations. 
#' #please keep only SNP mutations. This formulate just count the number, regardless of mutation type.
#'
#' @param MAF MAF files 
#' @param patients_col_name the column name of patients IDs, default is /strong{Tumor_Sample_Barcode}.
#' @param maxNum patients with mutations more than this number will be removed.
#' @param minNum patients with mutations less than this number will be removed.
#' @param patientsID giving the patients IDs and remaining mutations in these patients.
#' @return MAF after removing patients belongs to hyper-mutations 
#' @export

removing_patients= function(MAF, patients_col_name= "Tumor_Sample_Barcode",  maxNum= 3000, minNum=0,  patientsID = NULL){
  
  patients = data.frame( ID= rownames( table(MAF[,patients_col_name]) ), MutNum = as.numeric( table(MAF[,patients_col_name]) )  )
  
  message("Summary for Patients")
  
  print( summary(patients$MutNum) )
  
  message( sprintf("[1] the Total Patients are %s", nrow(patients) ))
  
  patients.num.left = patients[patients$MutNum >= minNum & patients$MutNum <= maxNum,  ]
  message( sprintf("[2] Patients remain %s  after mutations  cutoff between %s - %s ", nrow(patients.num.left) , minNum, maxNum))
  
  if(! is.null(patientsID)){
  patients.num.left = patients.num.left [ patients.num.left$ID %in% patientsID , ]  
  message( sprintf("[3] Patients remain %s  after Patients_ID cutoff ", nrow(patients.num.left) , maxNum, minNum))
  }
  
  message( sprintf("[4] %s patients are removing out", nrow(patients) - nrow(patients.num.left)  )   )
  
  print( summary(patients.num.left$MutNum) )
  
  mutation.num= nrow(MAF)
  
  MAF = MAF [ MAF[,patients_col_name] %in% patients.num.left$ID , ]
  
  message( sprintf("[5] In summary, %s mutations are removing out", mutation.num - nrow(MAF) )   )
  
  return(MAF)
}

