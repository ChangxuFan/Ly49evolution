rm(list = ls())
source("objects.R")
dir.create(work.dir, showWarnings = F, recursive = T)
samples <- SAMPLES.vivier.2022[9:10] # spleen samples only
sol <- utilsFanc::safelapply(samples, function(sample) {
  readRDS(paste0("vivier_2022/fast_check/sth/rna/per_sample/", sample, "/soi.Rds"))
}, threads = length(samples))
names(sol) <- samples

work.dir <- "vivier_2022/sth/merge_v2_sp/"
dir.create(work.dir, showWarnings = F, recursive = T)

soi <- cluster.pipe(soi = soi, assay = "RNA", is.sct.int = F,
                    pc.dim = 1:12, cluster.resolution = 0.5,
                    work.dir = work.dir, plot.dir = paste0(work.dir, "/plots/"), 
                    split.by = "sample",
                    project.name = "merge_v2_sp", save.rds = T, plot.common.markers = T,
                    do.sct = F)
