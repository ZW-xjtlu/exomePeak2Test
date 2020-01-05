#'@title Generate precision recall curve and related statistics directly from bam files in mayer's method.
#'@param bam_ip a vector containing the BAM files for IP samples.
#'@param bam_input a vector containing the BAM files for input samples.
#'@param txdb a TxDb object for the transcript annotation.
#'@param paired_end a logical value indicating wheather the data is paired end sequencing; default = FALSE.
#'@param bsgenome a BSgenome object for the genome sequence.
#'@param ground_truce_gr a GRanges for the ground truce of positive methylation sites, recommended to be in single based resolution.
#'@param N number of points sampled for each PRC curve; default = 200.
#'@param exp_label a character for the label of the experiment; default = "MeRIP_experiment_1".
#'@return The table for AUPRC and the AUROC curve will be saved on the disc under a folder named by exp_label.
#'
#'@examples
#' library(exomePeak2Test)
#' library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' library(exomePeak2)
#' library(BSgenome.Hsapiens.UCSC.hg19)
#'
#' whistle_predict <- readRDS("predictedResults.rds")
#' whistle_gr <- readRDS("exon_DRACH.rds")
#' whistle_gr$prob <- whistle_predict
#'
#' OneStep_PRC_Mayer(bam_ip = c("bam_local/SRR1035214_sorted.bam",
#'                        "bam_local/SRR1035216_sorted.bam",
#'                        "bam_local/SRR1035222_sorted.bam",
#'                        "bam_local/SRR1035224_sorted.bam"),
#'             bam_input = c("bam_local/SRR1035213_sorted.bam",
#'                           "bam_local/SRR1035215_sorted.bam",
#'                           "bam_local/SRR1035221_sorted.bam",
#'                           "bam_local/SRR1035223_sorted.bam"),
#'             txdb = TxDb.Hsapiens.UCSC.hg19.knownGene,
#'             bsgenome = Hsapiens,
#'             paired_end = FALSE,
#'             ground_truce_gr = whistle_gr[whistle_gr$prob > 0.5],
#'             exp_label = "MeRIP_experiment_1",
#'             N = 200)
#'
#'@import SummarizedExperiment GenomicRanges ggplot2
#'@export

OneStep_PRC_Mayer <- function(bam_ip,
                              bam_input,
                              txdb,
                              bsgenome = NULL,
                              paired_end = FALSE,
                              ground_truce_gr,
                              exp_label = "MeRIP_experiment_1",
                              N = 200){
  if(!dir.exists(exp_label)){
    dir.create(exp_label)
  }
  #saveRDS(sep, file.path(exp_label, paste0(file_name, ".rds"))) #latter use
  save_predict_peakcalling_mayer(BAM_IP = bam_ip,
                           BAM_INPUT = bam_input,
                           TXDB = txdb,
                           BSGNM = bsgenome,
                           PAIRED = paired_end,
                           TITLE = file.path(exp_label,"Temp"))

  result_lst = readRDS(file.path(exp_label,"Temp_mayer.rds"))
  gt_gr = ground_truce_gr
  group_name = exp_label

  pc_gr <- rowRanges( result_lst$dds )
  gt_gr <- subsetByOverlaps(gt_gr, pc_gr)
  pc_gr <- resize(pc_gr,201,fix = "center")

  pred_cfdr <- result_lst$cons_fdr_M
  pred_cfdr[is.na(pred_cfdr)] <- 0
  pred_cfdr[pred_cfdr == Inf] <- max(pred_cfdr[pred_cfdr != Inf])

  vec_pcc <- rep(NA,N)
  vec_recall <- rep(NA,N)

  cutoff <- quantile(pred_cfdr, seq(0,1,length.out = N), na.rm = TRUE)

  for(i in 1:N){
    pred_P <- pc_gr[rowSums( cbind( pred_cfdr >= cutoff[i] ) ) >= qbinom(0.95,ncol(cbind(pred_cfdr)),0.8)]
    vec_pcc[i] <- mean(pred_P %over% gt_gr)
    vec_recall[i] <- mean(gt_gr %over% pred_P)
  }

  vec_pcc <- c(vec_pcc, 1)
  vec_recall <- c(vec_recall, 0)

  plot_df <-  data.frame(recall = vec_recall,
                                       pcc = vec_pcc,
                                       value = "mayer",
                                       group = group_name )


  unlink(file.path(exp_label,"Temp_mayer.rds"))

  values <- unique( plot_df$value )
  groups <- unique( plot_df$group )

  #Fix no low points issue in combinded cutoff

    plot_df[1,"recall"] <- 1
    plot_df[1,"pcc"] <- min(  plot_df$pcc , na.rm = T)

  plot_df <- na.omit(plot_df)


  AUC <- function(x,y){
    require(zoo)
    id <- order(x)
    AUC <- sum(diff(x[id])*rollmean(y[id],2))
    return(AUC)
  }

  AUC_M <- matrix(NA,nrow = length(values), ncol = length(groups))
  rownames(AUC_M) = values
  colnames(AUC_M) = groups
  for(i in values){
    for(j in groups){
      df_sub <- plot_df[plot_df$value == i & plot_df$group == j,]
      AUC_M[i,j] <- AUC(df_sub$recall, df_sub$pcc)
    }
  }
  rm(df_sub)

  write.csv(AUC_M, file.path(exp_label, "AUC_M.csv"))

  write.csv(plot_df, file.path(exp_label, "plot_df.csv"))

}

get_pred_list_mayer <- function(SE_Peak_counts,
                                bsgenome,
                                txdb) {
  require(exomePeak2)
  require(DESeq2)
  pc_count_cutoff = 5
  bg_count_cutoff = 50
  binding_length = 25
  step_length = binding_length

  exome_bins_grl <- exomePeak2:::exome_bins_from_txdb(
    txdb = txdb,
    window_size = binding_length,
    step_size = step_length
  )

  if (is.null(bsgenome)) {} else {
    flanked_gr <- unlist(rowRanges(SE_Peak_counts))
    names(flanked_gr) <- gsub("\\..*$", "", names(flanked_gr))
    GC_freq <- as.vector(letterFrequency(Views(bsgenome, flanked_gr), letters = "CG"))
    sum_freq <- tapply(GC_freq, names(flanked_gr), sum)
    sum_freq <- sum_freq[names(rowRanges(SE_Peak_counts))]
    rowData(SE_Peak_counts)$gc_contents <-
      sum_freq / sum(width(rowRanges(SE_Peak_counts)))
    rm(flanked_gr, GC_freq, sum_freq)
  }

  #Bin width calculation
  rowData(SE_Peak_counts)$region_widths <-
    sum(width(rowRanges(SE_Peak_counts)))

  #Define background using prior knowledge of m6A topology


  rowData(SE_Peak_counts)$indx_bg = T

  #Change bins into initial widths
  rowData_tmp <- rowData(SE_Peak_counts)

  rowRanges(SE_Peak_counts) <- exome_bins_grl[as.numeric(rownames(SE_Peak_counts))]

  rowData(SE_Peak_counts) <- rowData_tmp

  rm(exome_bins_grl,rowData_tmp)

  SE_bins = SE_Peak_counts

  rm(SE_Peak_counts)

  count_cutoff = pc_count_cutoff

  design_IP_temp <- rep("input", ncol(SE_bins))

  design_IP_temp[colData(SE_bins)$design_IP] <- "IP"

  SE_bins$design_IP <- factor(design_IP_temp)

  SE_bins$design_IP <- relevel(SE_bins$design_IP, "input")

  rm(design_IP_temp)

  indx_count <- which(rowMeans(assay(SE_bins)) > count_cutoff)

  dds = DESeqDataSet(se = SE_bins[indx_count, ],
                     design = ~ design_IP)

  ######################################################
  #               Size factor estimation               #
  ######################################################
  off_sets = NULL

  if (!is.null(rowData(SE_bins)$gc_contents)) {
    dds$sizeFactor = 1

    indx_IP <- dds$design_IP == "IP"

    message("Estimating GC content correction factors for IP samples...")

    cqnObject_IP <-
      suppressMessages(
        cqn(
          assay(dds)[, indx_IP],
          lengths = rowData(dds)$region_widths,
          lengthMethod = "fixed",
          x = rowData(dds)$gc_contents,
          sizeFactors = dds$sizeFactor[indx_IP],
          verbose = FALSE,
          sqn = FALSE
        )
      )


    message("Estimating GC content correction factors for input samples...")

    cqnObject_input <-
      suppressMessages(
        cqn(
          assay(dds)[, !indx_IP],
          lengths = rowData(dds)$region_widths,
          lengthMethod = "fixed",
          x = rowData(dds)$gc_contents,
          sizeFactors = dds$sizeFactor[!indx_IP],
          verbose = FALSE,
          sqn = FALSE
        )
      )



    off_sets <- matrix(NA, nrow = nrow(dds), ncol = ncol(dds))

    rownames(off_sets) = rownames(dds)

    off_sets[, indx_IP] <- cqnObject_IP$offset

    off_sets[, !indx_IP] <- cqnObject_input$offset

    rm(cqnObject_IP, cqnObject_input, indx_IP)
  }


  N = sum(!dds$design_Treatment & dds$design_IP == "IP") * sum(!dds$design_Treatment & dds$design_IP == "input")

  fdr_M = matrix(NA,ncol = 1,nrow = nrow(dds))
  if(is.null(off_sets)){
  for (i in which(!dds$design_Treatment & dds$design_IP == "IP")){
    for( j in which(!dds$design_Treatment & dds$design_IP == "input")){
      fdr_M <-  cbind(fdr_M,mayers_hypertest(assay(dds)[,i],
                                             assay(dds)[,j]))
    }
  }
  }else{
    for (i in which(!dds$design_Treatment & dds$design_IP == "IP")){
      for( j in which(!dds$design_Treatment & dds$design_IP == "input")){
        fdr_M <-  cbind(fdr_M,mayers_hypertest(assay(dds)[,i],
                                               assay(dds)[,j],
                                               off_sets[,i],
                                               off_sets[,j]))
      }
    }
  }

  fdr_M <- fdr_M[,-1]

  return(list(
    dds = dds,
    cons_fdr_M = -log(fdr_M) #PS!
  ))

}

#' @rdname save_predict_peakcalling
#' @export

save_predict_peakcalling_mayer <- function(BAM_IP,BAM_INPUT,TXDB,BSGNM,PAIRED = FALSE,TITLE = "m6a"){

  require(magrittr)

  ReadCount <- get_reads_count(bam_ip = BAM_IP,
                               bam_input = BAM_INPUT,
                               txdb = TXDB,
                               paired_end = PAIRED)

  get_pred_list_mayer(SE_Peak_counts = ReadCount,
                      bsgenome = BSGNM,
                      txdb = TXDB)  %>% saveRDS(.,paste0(TITLE,"_mayer.rds"))

 }
