rm(list = ls())
source("objects.R")
library(DiffBind)
cage.signal.se.all <- readRDS("verify/B6/pass2/cage.b6.se.Rds")

### ATAC signal processing
# generate a summerizedExperiment (se) object with ATAC-seq signals at 
# peaks overlapping nanoCAGE-defined promoters.
# this part requires the AIAP pipeline output files. 
# So you need to run AIAP first...
# but I saved the resulting R object in a Rds file, so you can skip
# the ATAC signal processing part.
atac.insert.sites <- Sys.glob("path/to/my/aiap/pipeline/output/Processed*/step4.2*bedGraph")
names(atac.insert.sites) <- stringr::str_extract(string = atac.insert.sites, pattern = "\\d-D[np]") %>% 
  sub("(\\d)\\-(.+)", "\\2\\1", .)

atac.sample.df <- read.table("~/4dn/nk/fanc/atac2/diffbind/target_samples_mapped.tsv")[, c(1, 5)] %>% 
  `colnames<-`(c("SampleID", "Peaks"))

t.dba <- dba(minOverlap = 1, sampleSheet = atac.sample.df, scoreCol = 9, bRemoveM = T, 
             bRemoveRandom = T)

atac.peakset <- DiffBind::dba.peakset(DBA = t.dba, consensus = T, 
                                      bRemoveM = T, bRemoveRandom = T, minOverlap = 1,
                                      bRetrieve = T)

mcols( atac.peakset) <- NULL


atac.signal.df <- mclapply(names(atac.insert.sites), function(x) {
  file <- atac.insert.sites[x]
  signal <- gr.read.bedgraph(gr = atac.peakset, bdg = file)$score
  return(signal)
}, mc.cores = 6, mc.cleanup = T) %>% Reduce(cbind,.) %>% 
  as.data.frame() %>% `colnames<-`(names(atac.insert.sites))

# couldn't find where atac.signal.mat is from in peak_cor_2021-03-30.R. So I did this:
atac.signal.mat <- as.matrix(atac.signal.df)

atac.signal.se <- SummarizedExperiment(assays = list(counts = atac.signal.mat,
                                                     cpm = atac.signal.mat %*% diag(1/colSums(atac.signal.mat, na.rm = T)) * 1000000), 
                                       rowRanges = atac.peakset)

assays(atac.signal.se)$counts[is.na(assays(atac.signal.se)$counts)] <- 0

atac.signal.se <- se.fast.normalize(se = atac.signal.se, factor = 1000000, slot.name = "cpm")
saveRDS(atac.signal.se, "peak_cor/re_run/atac.signal.se.Rds")

### END: ATAC signal processing
atac.signal.se <- readRDS("atac.signal.se.Rds")

# generate the correspondence between ATAC peaks and nanoCAGE peaks. 
sync.o.all <- sync.atac.2.cage.2(cage.se = cage.signal.se.all, atac.se = atac.signal.se, genome = "mm10",
                                 bw.sample.map = "bw.sample.map.tsv", bw.ext.size =  500,
                                 threads = 6, rmbl = T, factor = 100000)
saveRDS(sync.o.all, "sync.o.all.Rds")
