#' Generate_RefCDS_sum
#'
#' Generate the RefCDS for the Ki/Ks model.
#'
#' @param RefCDS RefCDS
#' @param out saved file names
#'
#' @import rtracklayer
#' @import seqinr
#' @details load data before running this code.
#' @export

Generate_RefCDS_sum = function(RefCDS= RefCDS, out = out){

time_elapse <- function(start_time) {
  start_time <- as.POSIXct(start_time)
  dt <- difftime(Sys.time(), start_time, units="secs")
  # Since you only want the H:M:S, we can ignore the date...
  # but you have to be careful about time-zone issues
  format(.POSIXct(dt,tz="GMT"), "%H:%M:%S")
}

MutTable=list()

MutTable$N_EI = array(0, dim = c(192, 75))
MutTable$L_EI = array(0, dim = c(192, 75))

# Expanding the reference sequences [for faster access]
for (j in 1:length(RefCDS)) {
  RefCDS[[j]]$seq_cds = base::strsplit(as.character(RefCDS[[j]]$seq_cds), split="")[[1]]
  RefCDS[[j]]$seq_cds1up = base::strsplit(as.character(RefCDS[[j]]$seq_cds1up), split="")[[1]]
  RefCDS[[j]]$seq_cds1down = base::strsplit(as.character(RefCDS[[j]]$seq_cds1down), split="")[[1]]
  if (!is.null(RefCDS[[j]]$seq_splice)) {
    RefCDS[[j]]$seq_splice = base::strsplit(as.character(RefCDS[[j]]$seq_splice), split="")[[1]]
    RefCDS[[j]]$seq_splice1up = base::strsplit(as.character(RefCDS[[j]]$seq_splice1up), split="")[[1]]
    RefCDS[[j]]$seq_splice1down = base::strsplit(as.character(RefCDS[[j]]$seq_splice1down), split="")[[1]]
  }
}

############################################################
####### In this Part, fill MutTable$L_EI.
nt = c("A","C","G","T")
MutTable$L_EI = array(0, dim = c(192, 75))
time<-Sys.time()
for(j in 1:length(RefCDS)) {
  message(sprintf("run is %s", j))
  RefCDS[[j]]$L_EI=array(0, dim = c(192, 75))

  for(pos_ind in 1:RefCDS[[j]]$CDS_length ){
    for(MutNt in nt[-which(nt==RefCDS[[j]]$seq_cds[pos_ind])] ){
      ref3_cod = sprintf("%s%s%s", RefCDS[[j]]$seq_cds1up[pos_ind], RefCDS[[j]]$seq_cds[pos_ind], RefCDS[[j]]$seq_cds1down[pos_ind])
      mut3_cod = sprintf("%s%s%s", RefCDS[[j]]$seq_cds1up[pos_ind], MutNt, RefCDS[[j]]$seq_cds1down[pos_ind])
      codon_pos = c(ceiling(pos_ind/3)*3-2, ceiling(pos_ind/3)*3-1, ceiling(pos_ind/3)*3)
      old_codon = as.character(as.vector(RefCDS[[j]]$seq_cds[codon_pos]))
      pos_in_codon = pos_ind-(ceiling(pos_ind/3)-1)*3
      new_codon = old_codon; new_codon[pos_in_codon] = MutNt
      old_aa = seqinr::translate(old_codon)
      new_aa = seqinr::translate(new_codon)
      trisub = trinucsubsind[ paste(ref3_cod, mut3_cod, sep=">") ]
      aaid= which(AminoAcidDistance$AA== ifelse(new_aa<old_aa, sprintf("%s%s",new_aa, old_aa), sprintf("%s%s",old_aa, new_aa)) )
      MutTable$L_EI[trisub,aaid] = MutTable$L_EI[trisub,aaid] + 1 # Adding the mutation to the N matrices
      RefCDS[[j]]$L_EI[trisub,aaid] = RefCDS[[j]]$L_EI[trisub,aaid] + 1 # Adding the mutation to the N matrices

    }
  }
  message(sprintf("Time is %s", time_elapse(time)))
}

save(RefCDS,gr_genes, MutTable, file=sprintf("%s.rda",out) )

}
