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

  save_name <- gsub(".*\\/", "", pkc_rds_file)

  save_name <- gsub(".rds", "", save_name)

  code_readsGC <- call("plotReadsGC",
                         sep = expr_readRDS,
                         save_pdf_prefix = save_name
  ) %>% deparse #ellipse is not implemented yet, since call function cannot recognize it.

  index_last <- length(code_readsGC)

  code_readsGC[index_last] <- gsub(")$", "", code_readsGC[index_last] )

  code_readsGC[index_last] <- paste0(code_readsGC[index_last],
                                ",bsgenome = Hsapiens, txdb = TxDb.Hsapiens.UCSC.hg19.knownGene)")

  code_all <- c(code_library, code_readsGC)

  return(code_all)
}
