#' @title submit a bash commands on server background
#' @param command character vector of the R command used.
#' @param bashname the prefix name of files saved under bash.
#'
Rnohup <-
  function(command, bashname = "X") {

    writeLines(command,bashname)

    system(paste0("chmod +x ", bashname))

    outname <- gsub(".sh$",".out",bashname)

    system(paste0("nohup ", bashname, " > ", outname, " &"))

}
