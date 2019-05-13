#'@title Run peak callings in parallel on bash sever
#'@param coldata a \code{data.frame} which it columns containing the design information.
#'@param bam_dir the directory containing the bam files.
#'@param save_dir the directory to save the resulting RDS files of the \link{\code{summarizedExomePeaks}} object.
#'@param autocreate a logical, TRUE if create the \code{save_dir} when not finded.
#'@param parallel_num a numeric indicate the number of parallel process created; default = 6.
#'@param pc_function a peak calling function that takes the segragated coldata as input, and take a code string as output.
#'@param front_name a character string define the prefixes of the saved files.
#'@param ... arguments that passed to the peak calling function.
#'@details The coldata will be splitted only through the column of "Experiment".
#'@return multiple \link{\code{summarizedExomePeaks}} objects or others organized by experimental design.
#'@examples
#'ColData_hg19 <- read.csv("coldata_hg19.csv")
#'parallel_peak_calling(
#'coldata = ColData_hg19[ColData_hg19$Perturbation == "C",],
#'bam_dir = "~/merip_bam",
#'save_dir = "~/exomepeak2_test",
#'parallel_num = 9,
#'pc_function = exomePeak2_PC
#')
#'@import magrittr
#'@export
parallel_peak_calling <- function(coldata,
                                  bam_dir = "./",
                                  save_dir = "./",
                                  autocreate = T,
                                  parallel_num = 9,
                                  pc_function = exomePeak2_hg19,
                                  front_name = "pkc_",
                                  single_base = F,
                                  GC_correct = T,
                                  background = "Gaussian_mixture") {

  #check directories
  bam_dir <- gsub("/$", "", bam_dir)
  save_dir <- gsub("/$", "", save_dir)
  stopifnot(all(file.exists(
    paste0(bam_dir, "/", coldata$SRR_RUN, ".bam")
  )))
  if (!dir.exists(save_dir)) {
    stopifnot(autocreate)
    dir.create(save_dir)
  }

  #split coldata by experiments
  coldata_lst <- split(coldata, coldata$Experiment)

  #check the completeness of coldata
  index_complete <-
    sapply(coldata_lst, function(x)
      any(x$IP_input == "IP") & any(x$IP_input == "input"))

  if (any(!index_complete)) {
    warning(
      paste0(
        "The dataset ",
        paste(names(index_complete)[!index_complete], collapse = ", "),
        " have either no IP or input samples, the incomplete samples will not be analyzed."
      ),
      call. = FALSE,
      immediate. = TRUE
    )

    coldata_lst <- coldata_lst

  }

  #create codes with coding function
  code_lst <-
    lapply(coldata_lst,
           pc_function,
           bam_dir = bam_dir,
           front_name = front_name,
           single_base = single_base,
           GC_correct = GC_correct)

  #save R scripts on the system
  Rscript_names <- paste0(front_name, names(code_lst), ".R")

  mapply(function(x, y)
    writeLines(x, y),
    code_lst,
    paste0(save_dir, "/", Rscript_names))

  parallel_num = min(length(code_lst), parallel_num)

  #arrage the chunk of qsub according to parallel_num
  Rscript_commands <- paste0("Rscript ", save_dir, "/", Rscript_names)
  bash_chunks <-
    split(Rscript_commands, cut(seq_along(Rscript_commands),  parallel_num))
  bash_names <- paste0(front_name, seq_len(parallel_num), ".sh")
  names(bash_chunks) <- bash_names

  #Submit the bash commands with desired parallel number
  mapply(Rnohup,
         bash_chunks,
         paste0(save_dir, "/", bash_names))
}
