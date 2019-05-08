#'@title Produce code for guitar plot
#'@import magrittr
#'@import exomePeak2
#'@export
guitar_plot <- function(pkc_rds_file,
                          front_name,
                          ...){
  code_library <- c("library(exomePeak2)",
                    "library(TxDb.Hsapiens.UCSC.hg19.knownGene)",
                    "library(BSgenome.Hsapiens.UCSC.hg19)")

  expr_readRDS <- call("readRDS",
                       file = pkc_rds_file)

  code_GC_norm <- call("gcNormalization",
                       sep = expr_readRDS) %>% deparse

  arguments_annotation <- c("Hsapiens", "TxDb.Hsapiens.UCSC.hg19.knownGene")

  names(arguments_annotation) <- c("bsgenome", "txdb")

  code_GC_norm <- introduce_arg(code_GC_norm, arguments_annotation)

  code_GC_norm <- paste0("SEP <- ",code_GC_norm)

  save_name <- gsub(".*\\/", "", pkc_rds_file)

  save_name <- gsub(".rds", "", save_name)

  save_name <- paste0(front_name, save_name)

  code_plot <- call("plotGuitar",
                       save_pdf_prefix = save_name
  ) %>% deparse #ellipse is not implemented yet, since call function cannot recognize it.

  arguments_plot <- c("SEP","TxDb.Hsapiens.UCSC.hg19.knownGene")

  names(arguments_plot) <- c("sep", "txdb")

  code_plot<- introduce_arg( code_plot , arguments_plot )

  code_all <- c(code_library,code_GC_norm, code_plot)

  return(code_all)

}
