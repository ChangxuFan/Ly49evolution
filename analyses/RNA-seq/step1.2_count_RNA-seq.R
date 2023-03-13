source("~/R_packages/common/base.R")

bams <- Sys.glob("star_aligned/*filtered.bam")
bam.df <- data.frame(bamfile = bams,
                     sample = basename(bams) %>% gsub("_Aligned.+$|_RNA.+$", "", .))

liteRnaSeqFanc::rna.se.pe.featureCounts(
  bam.info = bam.df,
  annot = "../../annotations/gene_annotations/gencode.vM24.annotation.klraps.gtf",
  count.out.name = "featureCounts.tsv",
  thread = 8, count.frac.overlap = 0.5, other.count.options = "-p")

