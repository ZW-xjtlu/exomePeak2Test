#'@title Produce code for Peak Calling with exomePeak2
#'@imoort exomePeak2
#'@export
exomePeak2_PC <- function(coldata,
                                  front_name,
                                  ...){
  #Create specific representation of those code.
  code_library <- c("library(exomePeak2)",
                    "library(TxDb.Hsapiens.UCSC.hg19.knownGene)")

  expr_scanbam <- call("scanMeripBAM",
    bam_files = paste0(coldata$SRR_RUN, ".bam" ),
    design_ip = coldata$IP_input == "IP" ,
    paired_end = all(coldata$Lib == "Paired")
  )

  code_epcalling <- call("exomePeakCalling",
    merip_bams = expr_scanbam
  ) %>% deparse #ellipse is not implemented yet, since call function cannot recognize it.

  index_last <- length(code_epcalling)

  code_epcalling[index_last] <- gsub(")$", "", code_epcalling[index_last] )
  code_epcalling[index_last] <- paste0(code_epcalling[index_last],", txdb = TxDb.Hsapiens.UCSC.hg19.knownGene)")
  code_epcalling[1] <- paste0("SEP <- ", code_epcalling[1])

  code_save <- paste0("saveRDS( SEP, '", paste0(pre_name, coldata$Experiment[1]), ".rds' )")
  code_all <- c(code_library, code_epcalling, code_save)

  return(code_all)
}
