#'@title Produce code for Peak Calling with exomePeak2
#'@import magrittr
#'@import exomePeak2
#'@export
exomePeak2_hg19 <- function(
                          coldata,
                          bam_dir,
                          front_name,
                          single_base = FALSE,
                          GC_correct = TRUE,
                          correct_GC_bg = FALSE,
                          background = "Gaussian_mixture",
                          qtnorm = TRUE
                          ) {

  #Create specific representation of those code.
  code_library <- c("library(exomePeak2)",
                    "library(TxDb.Hsapiens.UCSC.hg19.knownGene)",
                    "library(BSgenome.Hsapiens.UCSC.hg19)")

  code_directory <- c(
    paste0("dir.create('", front_name, coldata$Experiment[1], "')"),
    paste0("setwd('", front_name, coldata$Experiment[1], "')" )
  )

  if(!single_base) {
   code_annot <- "sb_gr = NULL"
  }else{
   code_annot <- 'sb_gr = readRDS(system.file("extdata", "m6A_hg19_annot.rds", package = "exomePeak2"))'
  }

  bam_files <- paste0( bam_dir, "/", coldata$SRR_RUN, ".bam" )

  code_ep2 <- call("exomePeak2",
                       bam_ip = bam_files[coldata$IP_input == "IP"],
                       bam_input = bam_files[coldata$IP_input == "input"],
                       paired_end = all(coldata$Lib == "Paired"),
                       correct_GC_bg = correct_GC_bg,
                       background = background,
                       qtnorm = qtnorm,
                       export_format = "CSV",
                       save_plot_analysis = T
  ) %>% deparse

  arguments_plot <- c( ifelse(GC_correct,"Hsapiens","NULL"),
                       "TxDb.Hsapiens.UCSC.hg19.knownGene",
                       "sb_gr" )

  names(arguments_plot) <- c("bsgenome", "txdb", "mod_annot")

  code_ep2 <- introduce_arg( code_ep2 , arguments_plot )

  return(c(code_library, code_directory, code_annot, code_ep2))

}
