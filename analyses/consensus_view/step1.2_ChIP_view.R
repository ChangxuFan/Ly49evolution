source("objects.R")
abao.e1map1.2 <- readRDS("figure_map18/abaos/e1map1_refine//abao.Rds")

chip.df <- read.table("data/Tbet_Runx3_df.tsv", header = T)
rownames(chip.df) <- chip.df$track.name

abao.e1map1.2 <- aba.add.consensus(abao.e1map1.2)
abao.e1map1.2 <- aba.add.track(abao = abao.e1map1.2, 
                               track.df = chip.df[c(1,3), ], 
                               bTrack.df.regionID.regex = T,
                               prepend.regionID = T, threads = 1, track.type = "bdg")

aba.write.track(abao = abao.e1map1.2,  consensus.name = "cons_N_0.5", track.name = "bdg",
                track.type = "bdg", track.order = "aln", write.consensus = F, normalize = T,
                scale.to = 1000, add.constant = 1.5,
                out.dir = "figure_map18/pileup/chip_v1/bdg/", push.2.igv = T)
