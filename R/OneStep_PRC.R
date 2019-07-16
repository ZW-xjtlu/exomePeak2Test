#'@title Generate precision recall curve and related statistics directly from bam files.
#'@param bam_ip a vector containing the BAM files for IP samples.
#'@param bam_input a vector containing the BAM files for input samples.
#'@param txdb a TxDb object for the transcript annotation.
#'@param paired_end a logical value indicating wheather the data is paired end sequencing; default = FALSE.
#'@param bsgenome a BSgenome object for the genome sequence.
#'@param ground_truce_gr a GRanges for the ground truce of positive methylation sites, recommended to be in single based resolution.
#'@param N number of points sampled for each PRC curve; default = 200.
#'@param exp_label a character for the label of the experiment; default = "MeRIP_experiment_1".
#'@param glm_type type of GLM fit when peak calling; default = "DESeq2".
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
#' OneStep_PRC(bam_ip = c("bam_local/SRR1035214_sorted.bam",
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

OneStep_PRC <- function(bam_ip,
                        bam_input,
                        txdb,
                        bsgenome = NULL,
                        paired_end = FALSE,
                        ground_truce_gr,
                        exp_label = "MeRIP_experiment_1",
                        N = 200,
                        glm_type = "DESeq2"){
if(!dir.exists(exp_label)){
  dir.create(exp_label)
}
#saveRDS(sep, file.path(exp_label, paste0(file_name, ".rds"))) #latter use
save_predict_peakcalling(BAM_IP = bam_ip,
                         BAM_INPUT = bam_input,
                         TXDB = txdb,
                         BSGNM = bsgenome,
                         PAIRED = paired_end,
                         TITLE = file.path(exp_label,"Temp"),
                         GLM_TYPE = glm_type)

prc_df_nogc <- PRC_df(readRDS(file.path(exp_label,"Temp_noGC.rds")), ground_truce_gr, N = N, group_name = "no GC")
unlink(file.path(exp_label,"Temp_noGC.rds"))
prc_df_mclustgc <- PRC_df(readRDS(file.path(exp_label,"Temp_mclustGC.rds")), ground_truce_gr, N = N, group_name = "mclust GC")
unlink(file.path(exp_label,"Temp_mclustGC.rds"))
prc_df_unifgc <- PRC_df(readRDS(file.path(exp_label,"Temp_uniformGC.rds")), ground_truce_gr, N = N, group_name = "uniform GC")
unlink(file.path(exp_label,"Temp_uniformGC.rds"))
prc_df_priorgc <- PRC_df(readRDS(file.path(exp_label,"Temp_priorGC.rds")), ground_truce_gr, N = N, group_name = "prior GC")
unlink(file.path(exp_label,"Temp_priorGC.rds"))
prc_df_mclustgcqnorm <- PRC_df(readRDS(file.path(exp_label,"Temp_mclustGC_qtnorm.rds")), ground_truce_gr, N = N, group_name = "mclust GC qtnorm")
unlink(file.path(exp_label,"Temp_mclustGC_qtnorm.rds"))
prc_df_unifgcqnorm <- PRC_df(readRDS(file.path(exp_label,"Temp_uniformGC_qtnorm.rds")), ground_truce_gr, N = N, group_name = "uniform GC qtnorm")
unlink(file.path(exp_label,"Temp_uniformGC_qtnorm.rds"))
prc_df_priorgcqnorm <- PRC_df(readRDS(file.path(exp_label,"Temp_priorGC_qtnorm.rds")), ground_truce_gr, N = N, group_name = "prior GC qtnorm")
unlink(file.path(exp_label,"Temp_priorGC_qtnorm.rds"))

plot_df <- rbind(prc_df_nogc,
                 prc_df_mclustgc,
                 prc_df_unifgc,
                 prc_df_priorgc,
                 prc_df_mclustgcqnorm,
                 prc_df_unifgcqnorm,
                 prc_df_priorgcqnorm)

rm(prc_df_nogc,
   prc_df_mclustgc,
   prc_df_unifgc,
   prc_df_priorgc,
   prc_df_mclustgcqnorm,
   prc_df_unifgcqnorm,
   prc_df_priorgcqnorm)

values <- unique( plot_df$value )
groups <- unique( plot_df$group )

#Fix no low points issue in combinded cutoff
for(i in groups){
  indx_i <- min(which( plot_df$value == "combinded_cutoff" & plot_df$group == i ))
  plot_df[indx_i,"recall"] <- 1
  plot_df[indx_i,"pcc"] <- min(  plot_df$pcc[plot_df$value == "deseq_pvalue" & plot_df$group == i] )
}

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

plot_df$value <- factor(plot_df$value, levels = c("deseq_pvalue",
                                                  "deseq_lfc",
                                                  "c-test_pvalue",
                                                  "c-test_lfc",
                                                  "combinded_cutoff",
                                                  "random"))

p1 <- ggplot(plot_df) + geom_path(aes(x = recall, y = pcc, colour = value)) +
  geom_point(aes(x = recall, y = pcc, colour = value),size = 0.2) +
  theme_classic() +
  scale_color_brewer(palette = "Spectral") + labs(x = "recall", y = "precision") + facet_grid(~group)

ggsave(file.path(exp_label,"prc1.pdf"), p1, width = 15, height = 3)

p2 <- ggplot(plot_df) + geom_path(aes(x = recall, y = pcc, colour = group)) +
  geom_point(aes(x = recall, y = pcc, colour = group),size = 0.2) +
  theme_classic() +
  scale_color_brewer(palette = "Spectral") + labs(x = "recall", y = "precision") + facet_grid(~value)

ggsave(file.path(exp_label,"prc2.pdf"), p2, width = 15, height = 3)

p3 <- ggplot(plot_df[plot_df$value == "deseq_pvalue",]) + geom_path(aes(x = recall, y = pcc, colour = group)) +
  geom_point(aes(x = recall, y = pcc, colour = group),size = 0.2) +
  theme_classic() +
  scale_color_brewer(palette = "Spectral") + labs(x = "recall", y = "precision", title = "Precision Recall Curves of DESeq2 p values") + ylim(0.44,0.9)

ggsave(file.path(exp_label,"prc3.pdf"), p3, width = 6, height = 4)
}
