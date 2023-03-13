source("objects.R")
work.dir <- "vivier_2022/fast_check/sth/rna/"
plot.dir <- "vivier_2022/fast_check/plots/"
count.dir <- "vivier_2022/sth/count/"
paste0("mkdir -p ", plot.dir, " ", work.dir) %>% system()
samples <- SAMPLES.vivier.2022

sample.info <- data.frame(
  sample = samples,
  dir = paste0(count.dir, "/", samples, "/outs/filtered_feature_bc_matrix"),
  organ = sub("_rep.+$", "", samples)
)

write.table(x = sample.info, "vivier_2022/fast_check/sample_info.tsv",
            sep = "\t", quote = F, row.names = F, col.names = T)

sol <- sc.qc.construct.so(sample.info = sample.info,
                          project.name = "fast_check",
                          mt.pattern = "mt-")

sol <- sc.qc.metadata.find.elbow.core(sol = sol, meta.term = "nFeature_RNA",
                                      override = T)

t <- sc.rank.qc(sol = sol, feature = "nFeature_RNA", ymax = 10000,
                plot.dir = plot.dir, project.name = "ori")

lower.elbow <- rep(list(1000), length(samples)) %>% `names<-`(samples)
lower.elbow$SG_rep1 <- 1250
sol <- sc.qc.metadata.find.elbow.core(sol = sol, meta.term = "nFeature_RNA",
                                      elbow.list = lower.elbow,
                                      override = T)
t <- sc.rank.qc(sol = sol, feature = "nFeature_RNA", ymax = 10000,
                plot.dir = plot.dir, project.name = "lower_end")

sol <- sc.qc.elbow.filter(sol = sol,metas = c("nFeature_RNA"),
                          project.name = "fast_check", take.lower = F)
t <- sc.rank.qc(sol = sol, feature = "nFeature_RNA", ymax = 10000,
                plot.dir = plot.dir, project.name = "lower_end_post")

sol <- sc.qc.metadata.find.elbow.core(sol = sol, meta.term = "nFeature_RNA",
                                      override = T)
t <- sc.rank.qc(sol = sol, feature = "nFeature_RNA", ymax = 10000,
                plot.dir = plot.dir, project.name = "higher_end")

higher.elbow <- rep(list(2500), length(samples)) %>% `names<-`(samples)
higher.elbow$SG_rep1 <- 3750

sol <- sc.qc.metadata.find.elbow.core(sol = sol, meta.term = "nFeature_RNA",
                                      elbow.list = higher.elbow,
                                      override = T)
t <- sc.rank.qc(sol = sol, feature = "nFeature_RNA", ymax = 10000,
                plot.dir = plot.dir, project.name = "higher_end")


sol <- sc.qc.metadata.find.elbow.core(sol = sol, meta.term = "percent.mt",
                                      override = T)
t <- sc.rank.qc(sol = sol, feature = "percent.mt", ymax = 100, plot.dir = plot.dir)


sol <- sc.qc.elbow.filter(sol = sol,metas = c("percent.mt", "nFeature_RNA"),
                          project.name = "fast_check")

t <- sc.rank.qc(sol = sol, feature = "nFeature_RNA", ymax = 10000, plot.dir = plot.dir, project.name = "post_filter")
t <- sc.rank.qc(sol = sol, feature = "percent.mt", ymax = 100, plot.dir = plot.dir, project.name = "post_filter")


saveRDS(sol, paste0(work.dir, "/sol_v1.Rds"), compress = F)

sol <- readRDS( paste0(work.dir, "/sol_v1.Rds"))
sol <- utilsFanc::safelapply(sol, function(so) {
  root.name <- so@meta.data$sample[1]
  utilsFanc::t.stat(root.name)
  so <- cluster.pipe(soi = so, assay = "RNA", pc.dim = 1:8, cluster.resolution = 0.4, 
                     work.dir = paste0(work.dir, "/per_sample/", root.name),
                     plot.dir = paste0(work.dir, "/per_sample/", root.name, "/plots/"), 
                     project.name = root.name, metas.include = "sample", 
                     save.rds = T, do.sct = F,
                     plot.common.markers = T)
  return(so)
}, threads = length(samples))
