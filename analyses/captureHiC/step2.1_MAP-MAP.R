rm(list = ls())
source("objects.R")
hico.list <- readRDS("publication/hico.list.v1.Rds")

names(hico.list) <- sub("30", "_30", names(hico.list))
# get normalization factors from .hic files
# these represent the sequence coverage at each bin.
normo <- norm.vec.gen.m(.hic.vec.list = .HIC.list, genome = "mm10", chr = "chr6", threads = 8,
                        resolution.vec = c(1000, 5000, 10000, 25000))
dir.create("publication/SIPG/coverage/")
saveRDS(normo, "publication/SIPG/coverage/normo_v1.Rds")
# a copy is provided

enhancer.summits <- readRDS("publication/enhancer.summits.Rds")
# ATAC-seq peak summits for Ly49 MAPs

enhancer.bins <- enhancer.summits %>% utilsFanc::gr.fit2bin(bin.size = 5000, expand = T)
utilsFanc::write.zip.fanc(enhancer.bins, "publication/SIPG/enhancer.bins.bed")

enhancer.gi <- gr.full.interaction(enhancer.bins)
regions(enhancer.gi) <- gr.use.ly49(gr = regions(enhancer.gi))
enhancer.gi$int.name <- paste0(anchors(enhancer.gi, "first")$forth, "_", 
                               anchors(enhancer.gi, "second")$forth)
saveRDS(enhancer.gi, "publication/enhancer.gi.Rds")

trash <- SIPG.pipe(fg.gi = enhancer.gi, hico.list = hico.list, normo.list = normo, fg.ext = 15000, 
                   fg.id.col = "int.name", distance.wobble = 5000, resolutions = 5000, 
                   mapqs = c("inter.hic", "inter_30.hic"), norms = c("VC"),
                   out.dir = "publication/SIPG/grid_ext_15000_include_pseudo/", 
                   color.range = NULL, center.piece = NULL, debug = F, 
                   plot.only = F, plot.mat.grid = T,
                   plot.bar = T, replot.collapse = F, max.quantile = 0.98,
                   stat.formula = NULL)
# the output file publication/SIPG/grid_ext_15000_include_pseudo/inter.hic..res_5000..VC/mat_grid/mean_all.png
# is the matrix of interaction heatmaps for each combination of 2 MAPs.