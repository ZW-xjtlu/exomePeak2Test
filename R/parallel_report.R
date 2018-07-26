#'@title Run reports on peak callings in parallel with bash sever.
#'@param pkc_dir the directory containing the Peak calling results in RDS format.
#'@param pkc_prefix a character of the prefix of the peakcalling result files; default "pkc_".
#'@param pkc_suffix a character of the suffix of the peakcalling result files; default ".rds".
#'@param autocreate a logical, TRUE if create the \code{save_dir} when it does not exist; default TRUE.
#'@param save_dir the directory to save the reported files.
#'@param parallel_num an integer indicates the number of parallel process created; default = 6.
#'@param report_function a report function that takes the peak calling result as the input, and take a code string as output.
#'@param front_name a character string define the prefixes of the saved files.
#'@param ... arguments that can be passed to the report function.
#'@return multiple report files on peak calling results, which will use the same middle name.
#'@examples
#'parallel_report(
#'  pkc_dir = "~/exomepeak2_test"
#'  pkc_prefix = "pkc_",
#'  pkc_suffix = ".rds",
#'  save_dir = "~/exomepeak2_GC_peaks",
#'  parallel_num = 9,
#'  report_function = reads_GC_plot
#')
#'@import magrittr
#'@export

parallel_report <- function(
  pkc_dir = "~/exomepeak2_test",
  pkc_prefix = "pkc_",
  pkc_suffix = ".rds",
  autocreate = TRUE,
  save_dir = "~/exomepeak2_GC_peaks",
  parallel_num = 9,
  report_function = reads_GC_plot,
  front_name = "GCreads_"
){

  #check directories
  pkc_dir <- gsub("/$", "", pkc_dir)
  save_dir <- gsub("/$", "", save_dir)

  file_pkc_results <- grep( paste0(pkc_suffix, "$"), list.files(pkc_dir), value = T )

  file_pkc_results <- grep( pkc_prefix, file_pkc_results, value = T )

  stopifnot(length(file_pkc_results) != 0)

  if (!dir.exists(save_dir)) {
    stopifnot(autocreate)
    dir.create(save_dir)
  }

  #create the pkc files separately
  file_pkc_results <- paste0( pkc_dir, "/", file_pkc_results )

  #create codes with coding function
  code_lst <-
    lapply(file_pkc_results,
           report_function,
           front_name = front_name)

  #save R scripts on the system ( The following part is the same for all of the parallel functions. )
  Rscript_names <- paste0(front_name, names(code_lst), ".R")

  mapply(function(x, y)
    writeLines(x, y),
    code_lst,
    paste0(save_dir, "/", Rscript_names))

  parallel_num = min(length(code_lst), parallel_num)

  #arrage the chunk of qsub according to parallel_num
  Rscript_commands <- paste0("Rscript ", Rscript_names)
  bash_chunks <-
    split(Rscript_commands, cut(seq_along(Rscript_commands),  parallel_num))
  bash_names <- paste0(front_name, seq_len(parallel_num), ".sh")
  names(bash_chunks) <- bash_names

  #Submit the bash commands with desired parallel number
  mapply(Rnohup,
         bash_chunks,
         paste0(save_dir, "/", bash_names))

}
