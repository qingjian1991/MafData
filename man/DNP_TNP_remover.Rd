% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DNP_TNP_remover.R
\name{DNP_TNP_remover}
\alias{DNP_TNP_remover}
\title{DNP remover}
\usage{
DNP_TNP_remover(MAF, delete_recur = F)
}
\arguments{
\item{MAF}{MAF file, with relevant fields consisten with the GDC MAF specification (\url{https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/}).}

\item{delete_recur}{T/F statement about whether to remove recurrent (02) tumors.}
}
\description{
Function that removes likely dinucleotide variants that are actually labeled as single nucleotide variants.
 Also removes nucleotides 2 positions apart, as these are likely not single nucleotide events
 (analysis reveals that read counts for nucleotides two positions apart are almost
 perfectly correlated, meaning it is likely sequencing error).
 If the tumor sample names have TCGA naming convention and are
 under a column header named Tumor_Sample_Barcode,
 you can also automatically remove recurrent tumors (and get a warning
 about tumor types other than primary)
}
