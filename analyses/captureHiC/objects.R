source("~/R_packages/common/bioc.R")

suppressMessages(library(GenomicInteractions))
suppressMessages(library(circlize))
suppressMessages(library(ComplexHeatmap))
suppressMessages(library(ggpubr))
suppressMessages(library(Matrix))
suppressMessages({
  R.files <- Sys.glob("~/R_packages/v4c/R/*.R")
  for (R.file in R.files) {
    source(R.file)
  }
})


SAMPLES <- c("NK1-Dneg", "NK1-Dpos", "NK2-Dneg", "NK2-Dpos")

INTER.hic <- paste0("juicer/", SAMPLES, "_inter.hic") %>% normalizePath()
INTER30.hic <- paste0("juicer/", SAMPLES, "_inter_30.hic") %>% normalizePath()
## note: hic files are large. Due to size, they are available at 
# https://wangftp.wustl.edu/~cfan/Ly49pub/large_files/hic/

.HIC.list <- list(inter.hic = INTER.hic,
                  inter_30.hic = INTER30.hic)

SAMPLE.COLDATA <- data.frame(Ly49D = c("neg", "pos", "neg", "pos"),
                             bio.rep = c("brep1", "brep1", "brep2", "brep2"))
rownames(SAMPLE.COLDATA) <- SAMPLES
SAMPLE.COLDATA.inter <- SAMPLE.COLDATA
rownames(SAMPLE.COLDATA.inter) <- paste0(rownames(SAMPLE.COLDATA.inter), "_inter.hic")

SAMPLE.COLDATA.inter30 <- SAMPLE.COLDATA
rownames(SAMPLE.COLDATA.inter30) <- paste0(rownames(SAMPLE.COLDATA.inter30), "_inter_30.hic")

SAMPLE.COLDATA.hic <- rbind(SAMPLE.COLDATA.inter, SAMPLE.COLDATA.inter30)
SAMPLE.COLDATA.hic$mapq <- stringr::str_extract(rownames(SAMPLE.COLDATA.hic), "inter.+hic$")

REGION.large.v1 <- "chr6:129361438-130542791"


KLRA4.bins <- list()
KLRA4.bins$res_5000 <- list()
KLRA4.bins$res_5000$promoter <- "chr6:130065001-130070000" %>% utilsFanc::loci.2.gr()
KLRA4.bins$res_5000$promoter$forth <- "Ly49d"


# EL4.CMP.INTER.hic <- c("juicer/EL4_inter.hic",
#                        Sys.glob("subsample/sth/juicer/*2M_75bp/aligned/inter.hic"))
# EL4.CMP.INTER_30.hic <- c("juicer/EL4_inter_30.hic",
#                           Sys.glob("subsample/sth/juicer/*2M_75bp/aligned/inter_30.hic"))
# .HIC.list.EL4 <- list(inter.hic = EL4.CMP.INTER.hic,
#                       inter_30.hic = EL4.CMP.INTER_30.hic)
