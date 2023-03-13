rm(list = ls())
source("objects.R")
library(Rsamtools)
samples <- SAMPLES.vivier.2022[9:10] # spleen samples
Klra.genes <- paste0("Klra", c(1, 3:10, "11-ps", "13-ps", "14-ps", 17))
sample.map <- data.frame(
  sample = samples,
  dir = paste0("vivier_2022/sth/count/", samples, "/")
)

# note: Files are available at 
# https://wangftp.wustl.edu/~cfan/Ly49pub/large_files/vivier_2022_count.
# They are not on github because of size

soi <- readRDS("vivier_2022/sth/merge_v2_sp/soi.Rds")
soi$celltype <- "NK"
soi$celltype[soi$seurat_clusters %in% c("5", "6", "7", "8")] <- "others"

### QC: how many reads have the MM tag?
loci <- list()
loci$Ly49 <- "chr6:129814596-130398832"
extract.bam(so = soi, sample.map = sample.map, 
            bam.subset.range = loci$Ly49, xf.filter = 25, filter.thread = 2,
            cluster.ident = "celltype", group.ident = "sample",
            threads.each = 4, npar = 2,
            out.dir = "vivier_2022/sth/Ly49_QC/v1/extract_bam/")
# execute the resulting file (vivier_2022/sth/Ly49_QC/v1/extract_bam//run.sh)
# or use the resulting files that I provide...
# in bash
bam.dir <- "vivier_2022/sth/Ly49_QC/v1/extract_bam/"
Ly49.bams <- paste0(bam.dir, rep(c("NK", "others"), each = 2), "..",rep(samples, 2),
                    "_xf25.bam")
file.exists(unique(Ly49.bams))
sc.bam.qc(bams = Ly49.bams, bam.sample.names = rep(samples, 2), genes = Klra.genes, 
          out.dir = "vivier_2022/sth/Ly49_QC/v1/")
# For the spleen samples, Klra14-ps seems to be the only one with 8% MM. other 
# genes all below 1%.
### END: QC: how many reads have the MM tag?

########### MM correction:
sot <- sc.MM.correct(so = soi, MM.df = "vivier_2022/sth/Ly49_QC/v1/v1_MM_cell.tsv")
# spot check a few to see if they are correct:
t <- read.table("vivier_2022/sth/Ly49_QC/v1/v1_MM_cell.tsv", header = T, comment.char = "")
t.cell <- "Spleen_rep1#TACGCTCAGACCGCCT-1"
sot@assays$RNA@counts["Klra14-ps", t.cell] # 0
soi@assays$RNA@counts["Klra14-ps", t.cell] # 1
t %>% filter(cell == t.cell, GN == "Klra14-ps")
#        sample                           cell        GN n.UMI n.MM frac.MM
# 1 Spleen_rep1 Spleen_rep1#TACGCTCAGACCGCCT-1 Klra14-ps     1    1       1
rm(sot)
rm(t)
rm(t.cell)
soi <- sc.MM.correct(so = soi, MM.df = "vivier_2022/sth/Ly49_QC/v1/v1_MM_cell.tsv")
dir.create("vivier_2022/sth/merge_v3_MM/")
saveRDS(soi, "vivier_2022/sth/merge_v3_MM/soi.Rds")
########### END: MM correction;
