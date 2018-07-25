#' @title submit a bash commands on server background
#' @param command character vector of the R command used.
#' @param bashname the prefix name of files saved under bash.
#'
Rnohup <-
  function(command,bashname = "X") {
    writeLines(command,paste0(bashname,".sh"))
    system(paste0("chmod +x ./", paste0(bashname,".sh")))
    system(paste0("nohup ./",paste0(bashname,".sh"), ifelse(is.null(bashname),"",paste0(" > ",bashname,".out"))," &"))
  }
