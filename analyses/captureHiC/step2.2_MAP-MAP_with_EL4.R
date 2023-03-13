rm(list = ls())
source("objects.R")
SAMPLES <- paste0(SAMPLES, "_subsample")
SAMPLES <- c("EL4", SAMPLES)
EL4.CMP.INTER.hic <- paste0("juicer/withEL4/", SAMPLES, "_inter.hic")
file.exists(EL4.CMP.INTER.hic)

EL4.CMP.INTER_30.hic <- paste0("juicer/withEL4/", SAMPLES, "_inter_30.hic")
file.exists(EL4.CMP.INTER_30.hic)

.HIC.list.EL4 <- list(inter.hic = EL4.CMP.INTER.hic,
                      inter_30.hic = EL4.CMP.INTER_30.hic)
# these files are at https://wangftp.wustl.edu/~cfan/Ly49pub/large_files/hic_withEL4/, due to size


hico.list.ctrl <- hico.list.gen(inter.hic = EL4.CMP.INTER.hic, inter30.hic = EL4.CMP.INTER_30.hic, 
                                region = REGION.large.v1,
                                colData = DataFrame(sample = SAMPLES),
                                save.rds = "publication/hico.list.ctrl.Rds")
# hico.list.ctrl.Rds is offered

normo.ctrl <- norm.vec.gen.m(.hic.vec.list = .HIC.list.EL4, genome = "mm10", chr = "chr6", threads = 8,
                             resolution.vec = c(1000, 5000, 10000, 25000))
saveRDS(normo.ctrl, "publication/normo.ctrl.Rds")

enhancer.gi <- readRDS("publication/enhancer.gi.Rds")

trash <- SIPG.pipe(fg.gi = enhancer.gi, hico.list = hico.list.ctrl,
                   normo.list = normo.ctrl, fg.ext = 15000, 
                   fg.id.col = "int.name", distance.wobble = 5000, resolutions = 5000, 
                   mapqs = c("inter.hic", "inter_30.hic"), norms = c("VC"),
                   out.dir = "publication/SIPG/EL4_grid_ext_15000/", 
                   color.range = NULL, center.piece = c(4, 4), debug = F)

# plots will be generated later (step2.3)