library(exomePeak2Test)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(exomePeak2)
library(BSgenome.Hsapiens.UCSC.hg19)

f1 = system.file("extdata", "IP1.bam", package="exomePeak2")
f2 = system.file("extdata", "IP2.bam", package="exomePeak2")
f3 = system.file("extdata", "IP3.bam", package="exomePeak2")
f4 = system.file("extdata", "IP4.bam", package="exomePeak2")
IP_BAM = c(f1,f2,f3,f4)
f1 = system.file("extdata", "Input1.bam", package="exomePeak2")
f2 = system.file("extdata", "Input2.bam", package="exomePeak2")
f3 = system.file("extdata", "Input3.bam", package="exomePeak2")
INPUT_BAM = c(f1,f2,f3)
f2 = system.file("extdata", "mod_annot.rds", package="exomePeak2")
MOD_ANNO_GRANGE <- readRDS(f2)
GENE_ANNO_GTF = system.file("extdata", "example.gtf", package="exomePeak2")

whistle_gr <- readRDS("whistle_gr.rds")

#Human embryo sterm cell control
OneStep_PRC_Mayer(bam_ip = IP_BAM,
                  bam_input = INPUT_BAM ,
                  txdb = makeTxDbFromGFF(GENE_ANNO_GTF),
                  paired_end = FALSE,
                  ground_truce_gr = whistle_gr[whistle_gr$prob > 0.5],
                  exp_label = "hESC_C",
                  N = 200)

OneStep_PRC_Mayer(bam_ip = IP_BAM,
                  bam_input = INPUT_BAM ,
                  txdb = makeTxDbFromGFF(GENE_ANNO_GTF),
                  bsgenome = Hsapiens,
                  paired_end = FALSE,
                  ground_truce_gr = whistle_gr[whistle_gr$prob > 0.5],
                  exp_label = "hESC_C",
                  N = 200)
