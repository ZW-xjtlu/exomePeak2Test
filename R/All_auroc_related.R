get_reads_count <- function(bam_ip,
                            bam_input,
                            txdb,
                            paired_end = FALSE){
  require(exomePeak2)
  require(GenomicAlignments)
  require(BiocParallel)

  fragment_length = 100
  binding_length = 25
  step_length = binding_length
  merip_bams <- scanMeripBAM(
    bam_ip = bam_ip,
    bam_input = bam_input,
    paired_end = paired_end
  )

  exome_bins_grl <- exome_bins_from_txdb(
    txdb = txdb,
    window_size = binding_length,
    step_size = step_length
  )

  register(SerialParam())
  register(MulticoreParam(workers = 1))
  register(SnowParam(workers = 1))
  yieldSize(merip_bams) = 5000000 #Control the yield size to inhibit memory overflow

  SE_Peak_counts <- suppressWarnings( summarizeOverlaps(
    features = exomePeak2:::split_by_name(
      exomePeak2:::flank_on_exons(
        grl = exome_bins_grl,
        flank_length = fragment_length - binding_length,
        txdb = txdb,
        index_flank = FALSE
      )
    ),
    reads = merip_bams,
    param = Parameter(merip_bams),
    mode = "Union",
    inter.feature = FALSE,
    preprocess.reads = ifelse((LibraryType(merip_bams) == "1st_strand"),
                              exomePeak2:::reads_five_POS_rev,
                              exomePeak2:::reads_five_POS),
    singleEnd = !any(asMates(merip_bams)),
    ignore.strand = (LibraryType(merip_bams) == "unstranded"),
    fragments = any(asMates(merip_bams))
  ) )

  colData(SE_Peak_counts) <- DataFrame(metadata(merip_bams))

  return(SE_Peak_counts)
}

get_pred_list <- function(SE_Peak_counts,
                          txdb,
                          bsgenome = NULL,
                          background = "Gaussian_mixture",
                          qtnorm = TRUE){
  require(exomePeak2)
  require(DESeq2)
  require(cqn)
  glm_type = "DESeq2"
  pc_count_cutoff = 5
  bg_count_cutoff = 50
  correct_GC_bg = TRUE
  binding_length = 25
  step_length = binding_length

  exome_bins_grl <- exome_bins_from_txdb(
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

  #Count cutoff index
  rowData(SE_Peak_counts)$indx_gc_est <-
    rowMeans(assay(SE_Peak_counts)) >= bg_count_cutoff

  #Model based clustering
  m6A_prior = F

  if (background == "Gaussian_mixture") {
    message("Identifying background with Gaussian Mixture Model...")

    rowData(SE_Peak_counts)$indx_bg <- exomePeak2:::mclust_bg(se_peak_counts = SE_Peak_counts)

    if (sum(rowData(SE_Peak_counts)$indx_gc_est &
            rowData(SE_Peak_counts)$indx_bg) < 30) {

      warning(
        "Background bin # < 30 using mclust, search background with m6A-seq prior",
        call. = FALSE,
        immediate. = FALSE
      )
      m6A_prior = T
    }

  }


  #Define background using prior knowledge of m6A topology

  if (background == "m6Aseq_prior" | m6A_prior) {
    indx_UTR5 <-
      rowRanges(SE_Peak_counts) %over% fiveUTRsByTranscript(txdb)

    indx_longexon <-
      rowRanges(SE_Peak_counts) %over% exons(txdb)[width(exons(txdb)) >= 400]

    rowData(SE_Peak_counts)$indx_bg = !(indx_UTR5 | indx_longexon)

    rm(indx_UTR5, indx_longexon, m6A_prior)
  }

  if (background == "all") {
    rowData(SE_Peak_counts)$indx_bg = T
  }

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

  dds$sizeFactor = estimateSizeFactorsForMatrix(assay(dds))

  if (!is.null(rowData(SE_bins)$gc_contents)) {
    indx_IP <- dds$design_IP == "IP"

    message("Estimating GC content correction factors for IP samples...")

    if(correct_GC_bg) {
      subindex = which(rowData(dds)$indx_bg & rowData(dds)$indx_gc_est)
    }else{
      subindex = which(rowData(dds)$indx_gc_est)
    }

    cqnObject_IP <-
      suppressMessages(
        cqn(
          assay(dds)[, indx_IP],
          lengths = rowData(dds)$region_widths,
          lengthMethod = "fixed",
          x = rowData(dds)$gc_contents,
          sizeFactors = dds$sizeFactor[indx_IP],
          subindex = subindex,
          verbose = FALSE,
          sqn = qtnorm
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
          subindex = subindex,
          verbose = FALSE,
          sqn = qtnorm
        )
      )


    rm(subindex)

    glm_off_sets <- matrix(NA, nrow = nrow(dds), ncol = ncol(dds))

    rownames(glm_off_sets) = rownames(dds)

    glm_off_sets[, indx_IP] <- cqnObject_IP$glm.offset

    glm_off_sets[, !indx_IP] <- cqnObject_input$glm.offset

    rm(cqnObject_IP, cqnObject_input, indx_IP)

    #normalization to make the row geometric means = 0 (since DESeq2 only cares about the difference).
    #and this norm factor is still under the original scale (not log scale glm off set).
    centered_off_sets <-
      exp(glm_off_sets) / exp(rowMeans(glm_off_sets))

    normalizationFactors(dds) <- centered_off_sets

    rm(glm_off_sets, centered_off_sets)
  }


  ######################################################
  #             Generalized Linear Model               #
  ######################################################

  if (glm_type == "Poisson") {
    dispersions(dds) = 0
  }

  if (glm_type == "NB") {
    dds = estimateDispersions(dds, fitType = "mean")
  }

  if (glm_type == "DESeq2") {
    dds = suppressMessages( estimateDispersions(dds) )
  }

  dds = suppressMessages( nbinomWaldTest(dds) )

  res <- results(dds, altHypothesis = "greater")

  if(is.null(normalizationFactors(dds))){
    sizeFactor_M <- matrix(dds$sizeFactor,nrow = nrow(dds), ncol = ncol(dds), byrow = T)
  }else{
    sizeFactor_M <- normalizationFactors(dds)
  }

  N = sum(!dds$design_Treatment & dds$design_IP == "IP") * sum(!dds$design_Treatment & dds$design_IP == "input")

  pvalue_M = matrix(NA,ncol = 1,nrow = nrow(dds))
  fdr_M = matrix(NA,ncol = 1,nrow = nrow(dds))
  logFC_M = matrix(NA,ncol = 1,nrow = nrow(dds))

  for (i in which(!dds$design_Treatment & dds$design_IP == "IP")){
    for( j in which(!dds$design_Treatment & dds$design_IP == "input")){
      test_result <- exomePeak2:::ctest(assay(dds)[,i],
                                        assay(dds)[,j],
                                        sizeFactor_M[,i],
                                        sizeFactor_M[,j],
                                        fold = 1)
      fdr_M <-  cbind(fdr_M,test_result$fdr)
      logFC_M <- cbind(logFC_M,test_result$log2FC)
      pvalue_M <- cbind(logFC_M,test_result$pvalue)
    }
  }

  logFC_M <- logFC_M[,-1]
  fdr_M <- fdr_M[,-1]

  return(list(
    dds = dds,
    DESeq2Results = res,
    cons_pvalue_M = pvalue_M,
    cons_log2FC_M = logFC_M,
    cons_fdr_M = fdr_M
  ))

}

#' @rdname save_predict_peakcalling
#' @export

save_predict_peakcalling <- function(BAM_IP,BAM_INPUT,TXDB,BSGNM,PAIRED = FALSE,TITLE = "m6a"){

  require(magrittr)

  ReadCount <- get_reads_count(bam_ip = BAM_IP,
                               bam_input = BAM_INPUT,
                               txdb = TXDB,
                               paired_end = PAIRED)

  get_pred_list(SE_Peak_counts = ReadCount,
                txdb = TXDB,
                bsgenome = NULL,
                background = "all",
                qtnorm = FALSE)  %>% saveRDS(.,paste0(TITLE,"_noGC.rds"))

  get_pred_list(SE_Peak_counts = ReadCount,
                txdb = TXDB,
                bsgenome = BSGNM,
                background = "Gaussian_mixture",
                qtnorm = FALSE)  %>% saveRDS(.,paste0(TITLE,"_mclustGC.rds"))

  get_pred_list(SE_Peak_counts = ReadCount,
                txdb = TXDB,
                bsgenome = BSGNM,
                background = "Gaussian_mixture",
                qtnorm = TRUE)  %>% saveRDS(.,paste0(TITLE,"_mclustGC_qtnorm.rds"))

  get_pred_list(SE_Peak_counts = ReadCount,
                txdb = TXDB,
                bsgenome = NULL,
                background = "m6Aseq_prior",
                qtnorm = FALSE)  %>% saveRDS(.,paste0(TITLE,"_priorGC.rds"))

  get_pred_list(SE_Peak_counts = ReadCount,
                txdb = TXDB,
                bsgenome = NULL,
                background = "m6Aseq_prior",
                qtnorm = TRUE)  %>% saveRDS(.,paste0(TITLE,"_priorGC_qtnorm.rds"))

  get_pred_list(SE_Peak_counts = ReadCount,
                txdb = TXDB,
                bsgenome = BSGNM,
                background = "all",
                qtnorm = FALSE)  %>% saveRDS(.,paste0(TITLE,"_uniformGC.rds"))

  get_pred_list(SE_Peak_counts = ReadCount,
                txdb = TXDB,
                bsgenome = BSGNM,
                background = "all",
                qtnorm = TRUE)  %>% saveRDS(.,paste0(TITLE,"_uniformGC_qtnorm.rds"))

}
