#' @title Introduce argument into a code chunk.
#' @param code a character vector containing the code.
#' @param arguments a character vector with names being the arguments and entries being the inputs.
#' @export
introduce_arg<- function(code, arguments) {
  index_last <- length(code)

  code[index_last] <- gsub(")$", "", code[index_last] )

  argument_text <- paste( paste0(names(arguments)," = ",arguments) , collapse = ", ")

  code[index_last] <- paste0(code[index_last],
                                    ", ",argument_text,")")

  return(code)
}
