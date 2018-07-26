#'@title Produce code for Peak Calling with exomePeak2
#'@import magrittr
#'@import exomePeak2
#'@export
reads_GC_plot <- function(pkc_rds_file,
                          front_name,
                          ...){
  #Create specific representation of those code.
  code_library <- c("library(exomePeak2)",
                    "library(TxDb.Hsapiens.UCSC.hg19.knownGene)",
                    "library(BSgenome.Hsapiens.UCSC.hg19)")

  expr_readRDS <- call("readRDS",
                       file = pkc_rds_file)

  code_GC_norm <- call("GCnormalization",
                       sep = expr_readRDS) %>% deparse

  arguments_annotation <- c("Hsapiens", "TxDb.Hsapiens.UCSC.hg19.knownGene")

  names(arguments_annotation) <- c("bsgenome", "txdb")

  code_GC_norm <- introduce_arg(code_GC_norm, arguments_annotation)

  save_name <- gsub(".*\\/", "", pkc_rds_file)

  save_name <- gsub(".rds", "", save_name)

  save_name <- paste0(front_name, save_name)

  code_readsGC <- call("plotReadsGC",
                         save_pdf_prefix = save_name
  ) %>% deparse #ellipse is not implemented yet, since call function cannot recognize it.

  arguments_plot <- code_GC_norm

  names(arguments_plot) <- "sep"

  code_readsGC <- introduce_arg( code_readsGC , arguments_plot )

  code_all <- c(code_library, code_readsGC)

  return(code_all)
}
