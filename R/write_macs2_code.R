#'@title Write code for macs2 (only for auprc purpose).
#'@param bam_ip a vector containing the BAM files for IP samples.
#'@param bam_input a vector containing the BAM files for input samples.
#'@param name a character for the name of the output sh file and peak calling file.
#'@return A file named by \code{name}.
#'@examples
#' OneStep_PRC(bam_ip = c("bam_local/SRR1035214_sorted.bam",
#'                        "bam_local/SRR1035216_sorted.bam",
#'                        "bam_local/SRR1035222_sorted.bam",
#'                        "bam_local/SRR1035224_sorted.bam"),
#'             bam_input = c("bam_local/SRR1035213_sorted.bam",
#'                           "bam_local/SRR1035215_sorted.bam",
#'                           "bam_local/SRR1035221_sorted.bam",
#'                           "bam_local/SRR1035223_sorted.bam"),
#'             name = "hNPC_control")
#'
#'@export

write_macs2_code <- function( bam_ip, bam_input, name ) {

 code_lines  <- c("macs2 callpeak -t",
                  bam_ip,
                  "-c",
                  bam_input,
                  paste0("-f BAM -g hs -n ", name, " -B -p 0.9999"))

 code_lines <- paste(code_lines, "\\")

 writeLines(code_lines, con = paste0(name, ".sh"))

 system(paste0("chmod u+x ",name, ".sh"))

}
