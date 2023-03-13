source("objects.R")
# a rough alignment that covers entire Ly49 gene bodies
aba.df <- readRDS("data/aba_df_full_v2.Rds")
abao.full <- aba.create(aba.df = aba.df, df.is.bed = F, work.dir = "figure_map18/abaos/full_v2/", 
                       threads = 8)

# subset the alignment to only contain exon 1 and upstream MAP1/8 regions, based on 
# B6.Ly49c.
map1.Ly49c <- abao.full@ori.df %>% 
  dplyr::filter(regionID == "B6.Ly49c")
map1.Ly49c$start <- 130337470
map1.Ly49c$end <- map1.Ly49c$start + 4974 # where synteny broke between a/g and c/i

abao.e1map1 <- aba.subset(abao = abao.full, 
                          smart.cut.df = map1.Ly49c, regionIDs.exclude = c("Ly49q\\d*"), 
                          new.work.dir = "figure_map18/abaos/e1map1/")
# filter out Ly49s that don't align well (mostly due to issues of the assemblies)
map1.df <- abao.e1map1@ori.df %>% 
  dplyr::filter(regionID != "129.Ly49u") %>% 
  dplyr::filter(!grepl("Ly49q\\d*$", regionID)) %>% 
  dplyr::filter(regionID != "NOD.Ly49p1")

map1.df[map1.df$regionID == "129.Ly49ui", "start"] <- 433431
abao.e1map1.2 <- aba.create(aba.df = map1.df, df.is.bed = F, 
                            work.dir = "figure_map18/abaos/e1map1_refine/", threads = 8)

# known abnormalities:
# gap: MAP1: 129.Ly49ec1; NOD.Ly49p1; 129.Ly49ui; MAP8: 129.Ly49r, NOD.Ly49p1
# no MAP1: Ly49x (B6 and NOD); 129.Ly49lr; 
# incomplete MAP1: NOD.Ly49p2
# make sure this object contains known Ly49 genes:
strains <- c("B6", "129", "NOD")
missing.Ly49s <- lapply(strains, function(strain) {
  if (strain == "B6")
    genes <- readLines(paste0("~/genomes/ly49/gene_names/", "mm10", ".txt"))
  else
    genes <- readLines(paste0("~/genomes/ly49/gene_names/", strain, ".txt"))
  genes <- paste0(strain, ".", genes)
  not.found <- genes %>% .[!.%in%abao.e1map1.2@ori.df$regionID]
  return(not.found)
})

names(missing.Ly49s) <- strains
missing.Ly49s

# subset out MAP8 to check alignment quality:
map8.df <- readRDS("data/aba_df_map8.Rds") %>% 
  dplyr::filter(regionID == "B6.Ly49c")
aba.subset("figure_map18/abaos/e1map1_refine//abao.Rds", smart.cut.df = map8.df,
           new.work.dir = "figure_map18/abaos/e1map1_sub_map8/")

# add and write annotation tracks:
abao.e1map1.2 <- aba.add.track(abao.e1map1.2, 
                               track.file = "data/anno_B6.Ly49c_simple.bed",
                               track.regionID.regex = "B6", track.name = "bed",
                               prepend.regionID = T, threads = 1, track.type = "bed")

aba.write.track(abao = abao.e1map1.2, consensus.name = "cons_N_0.5",
                write.consensus = T, consensus.add.N = 10000,
                fill.gap = T, broadcast = T, 
                track.name = "bed", track.type = "bed", 
                out.dir = "figure_map18/pileup/v1/bed/")

# add and write ATAC-seq signal tracks
track.df <- read.table("data/atac_track_df.tsv", header = T)
track.df.q8 <- track.df %>% 
  mutate(track.file = sub("_AS0", "_AS0_mapq8", track.file))
track.list <- list(noMAPQ = track.df, MAPQ8 = track.df.q8)
lapply(names(track.list), function(type) {
  df <- track.list[[type]]
  abao.e1map1.2 <- aba.add.track(abao = abao.e1map1.2, track.df = df,
                                 bTrack.df.regionID.regex = T,
                                 prepend.regionID = T, threads = 1, track.type = "bdg")
  
  aba.write.track(abao = abao.e1map1.2,  consensus.name = "cons_N_0.5", track.name = "bdg",
                  track.type = "bdg", track.order = "aln", write.consensus = F,
                  normalize = T, add.constant = 1.5,
                  out.dir = paste0("figure_map18/pileup/v1/atac_", type, "/"), push.2.igv = T)
})
