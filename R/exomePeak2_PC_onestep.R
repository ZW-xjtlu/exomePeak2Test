#'@title Produce code for Peak Calling with exomePeak2
#'@import magrittr
#'@import exomePeak2
#'@export
exomePeak2_PC_onestep <- function(coldata,
                          bam_dir,
                          front_name,
                          ...){
  #Create specific representation of those code.
  code_library <- c("library(exomePeak2)",
                    "library(TxDb.Hsapiens.UCSC.hg19.knownGene)",
                    "library(BSgenome.Hsapiens.UCSC.hg19)")

  code_directory <- c( paste0("dir.create('", front_name, coldata$Experiment[1], "')"),
                       paste0("setwd('", front_name, coldata$Experiment[1], "')" )
  )

  bam_files <- paste0( bam_dir, "/", coldata$SRR_RUN, ".bam" )

  code_ep2 <- call("exomePeak2",
                       bam_ip = bam_files[coldata$IP_input == "IP"],
                       bam_input = bam_files[coldata$IP_input == "input"],
                       paired_end = all(coldata$Lib == "Paired"),
                       export_format = "RDS"
  ) %>% deparse

  arguments_plot <- c("Hsapiens","TxDb.Hsapiens.UCSC.hg19.knownGene")

  names(arguments_plot) <- c("bsgenome", "txdb")

  code_ep2 <- introduce_arg( code_ep2 , arguments_plot )

  return(c(code_library, code_directory,  code_ep2))
}
