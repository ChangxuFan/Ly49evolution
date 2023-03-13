source("objects.R")
work.dir <- "de"
dir.create(work.dir)
rawmat <- liteRnaSeqFanc::featureCounts.2.mat(in.tsv = "featureCounts.tsv.renamed")
rownames(rawmat) <- utilsFanc::klra.ly49.translate(rownames(rawmat), bKeep.not.found = T)
colnames(rawmat) <- colnames(rawmat) %>% sub("X(.+)_(.+)$", "\\2_\\1", .)
colnames(rawmat)
# [1] "WT_601.2" "KO_601.3" "KO_601.4" "WT_601.5"

### these correspond to the following libraries from GEO:
# MAP8.B6.Ly49m_WT_M6_RNA_rep1; MAP8.B6.Ly49m_KO_M6_RNA_rep1;
# MAP8.B6.Ly49m_KO_M6_RNA_rep2; MAP8.B6.Ly49m_WT_M6_RNA_rep2
###

sample.info <- data.frame(sample = colnames(rawmat)) %>% 
  dplyr::mutate(genotype = sub("_.+$", "", sample),
                mouse = sub("^.+_", "", sample),
                celltype = "NK")
cutoff <- 5
n.pass <- 2
bMat <- rawmat > cutoff
sum(rowSums(bMat) >= n.pass)
# 12700.
filtered.mat <- rawmat[rowSums(bMat) >= n.pass,]
rowMins(filtered.mat) %>% summary()
rowMaxs(filtered.mat) %>% summary()

# > rowMaxs(filtered.mat) %>% summary()
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#       6      34     130     354     355   28021 
# > rowMins(filtered.mat) %>% summary()
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#       0      21      98     275     280   19129 

### note for reader: the following below was originally written for ATAC-seq,
# with the goal to get rid of DiffBind from my analyses pipeline. I ended up
# using it mostly for RNA-seq... But the name stayed.
# this is essentially a wrapper for DESeq2.
de <- diffbind.2.raw.a2bl(sample.info = sample.info,
                          mat = filtered.mat,
                          cluster.ident = "celltype",sample.col = "sample", 
                          filter.nz = F, filter.size = c(0, 0.99), sequential.filter = T,
                          coldata.columns = c("genotype", "mouse", "celltype"), 
                          design.formula = ~ genotype, 
                          contrast = c("genotype", "KO", "WT"), 
                          threads = 1,
                          work.dir = work.dir)

### make Ly49 expression bar plot
gtf <- rtracklayer::import("../../annotations/gene_annotations/Ly49_mm10.gtf")
gtf <- gtf[gtf$type == "gene"]
gtf$gene <- gtf$gene_name %>% sub(".*Ly49", "", .)
gtf <- resize(gtf, width = 1, fix = "center")
gtf.df <- gtf %>% as.data.frame() %>% dplyr::select(gene, start) %>% 
  dplyr::rename(coord = start) 

gtf.df <- gtf.df %>% dplyr::group_by(gene) %>% 
  dplyr::summarise(coord = floor(mean(coord))) %>% 
  dplyr::ungroup() %>% as.data.frame() %>% arrange(coord)
write.table(gtf.df, "Ly49_coord_simple.tsv", quote = F,
            sep = "\t", col.names = T, row.names = F)

gtf.df <- read.table("Ly49_coord_simple.tsv", header = T)

metadata(de$NK$res)$filterThreshold
#       0% 
# 3.468484
# acceptable
Ly49.expr <- paste0("Ly49", c("f", "d", "k", "h", "n", "i", "g", "j", "m", "c", "a"))
bar.df <- de.plot.bar(de = de, genes = Ly49.expr, clusters = "NK", return.df = T)  
q.df <- de$NK$res.exp %>% filter(gene %in% Ly49.expr) %>% 
  dplyr::select(padj, gene) %>% 
  dplyr::mutate(group.1 = "WT", group.2 = "KO", gene = sub("Ly49", "", gene)) %>% 
  dplyr::rename(p = padj)
dir.create(paste0(work.dir, "/plot_bar/"))
p <- utilsFanc::barplot.pub.3(df = bar.df, x = "gene", y = "expr", color.by = "genotype",
                   genomic.x = gtf.df, spread.bin.size = 100,
                   bar.width = 0.5, dodge.width = 0.6,
                   add.pval = T, pval.group.1 = "WT", pval.group.2 = "KO", 
                   external.p.df = q.df) %>% 
  utilsFanc::theme.fc.1() +
  ggsave(paste0(work.dir, "/plot_bar/Ly49_expressed_bar_2.pdf"), device = cairo_pdf, 
         width = 4.5, height = 1.0)

### END: make Ly49 expression bar plot