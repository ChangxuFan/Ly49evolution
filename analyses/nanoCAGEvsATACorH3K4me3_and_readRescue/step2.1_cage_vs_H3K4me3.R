rm(list = ls())
source("objects.R")
k4me3.sample.df <- data.frame(bamfile = c("/path/to/K4me3_NK_WT_exvivo_rep1_chip_GSM4314396.bam",
                                          "/path/to/K4m3_NK_WT_exvivo_rep2_chip_GSM4314407.bam"),
                              sample = c("rep1", "rep2"))

cage.b6.se <- readRDS("verify/B6/pass2/cage.b6.se.Rds")
k4me3.se <- create.histone.se(gr = rowRanges(cage.b6.se), buffer.up = 1000, buffer.down = 3000, 
                              samples.df = k4me3.sample.df,  bSingleEnd = T, 
                              out.dir = "cage_histone/k4me3/", root.name = "k4me3",
                              scale.factor = 10000)

k4me3.synco <- synco.gen.ez(cage.se = cage.b6.se, atac.se = k4me3.se, use.cage.as.j = T)

k4me3.synco <- synco.prep(synco = k4me3.synco, annotate = T, rmbl = T, genome = "mm10")
dir.create("cage_histone/k4me3_2023")
saveRDS(k4me3.synco, "cage_histone/k4me3_2023/k4me3.synco.Rds")
saveRDS(k4me3.se, "cage_histone/k4me3_2023/k4me3.se.Rds")
k4me3.cor <- synco.range.cor(synco = k4me3.synco, outfile = "cage_histone/k4me3_2023/k4me3_range_cor.pdf",
                             add.cor = T, use.template.2 = T)

