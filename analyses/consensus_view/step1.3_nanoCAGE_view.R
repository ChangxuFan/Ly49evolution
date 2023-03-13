source("objects.R")

# this is what I think are the high quality ones after manual inspection of alignment
Ly49.cage.hq.2 <- c("129.Ly49i1",
                    "B6.Ly49i",
                    "NOD.Ly49i1",
                    "B6.Ly49m",
                    "NOD.Ly49m",
                    "B6.Ly49h",
                    "NOD.Ly49h",
                    "NOD.Ly49u",
                    "NOD.Ly49d",
                    "129.Ly49r",
                    "B6.Ly49c",
                    "B6.Ly49a",
                    "129.Ly49p",
                    "129.Ly49o",
                    "NOD.Ly49a",
                    "B6.Ly49g",
                    "129.Ly49g",
                    "NOD.Ly49e",
                    "B6.Ly49f",
                    "129.Ly49lr",
                    "B6.Ly49d",
                    "NOD.Ly49g2",
                    "129.Ly49v")

# revise this list by expression level
expr.non.pd <- lapply(c("B6", "129", "NOD"), function(strain) {
  x <- ifelse(strain == "129", "m129", strain)
  expr.df <- read.table(paste0("data/nanoCAGE_expr_ranked/", x, ".tsv"), header = T)
  exp.col <- ifelse(strain == "B6", 4, 2)
  expr.df <- expr.df[expr.df[, exp.col] > 50 & !expr.df$is.pseudo,]
  genes <- paste0(strain, ".", expr.df$gene)
  return(genes)
}) %>% unlist()



expr.non.pd %>% .[!. %in% Ly49.cage.hq.2]
# "NOD.Ly49p3" "NOD.Ly49p1"
Ly49.cage.hq.2 %>% .[!. %in% expr.non.pd]
# [1] "B6.Ly49m"   "129.Ly49lr"
# decision: just add p1, p3 in and remove lr and m

# revised list
Ly49.cage.hq.3 <- c("129.Ly49i1",
                    "B6.Ly49i",
                    "129.Ly49r",
                    "B6.Ly49c",
                    "B6.Ly49h",
                    "NOD.Ly49d",
                    "NOD.Ly49h",
                    "NOD.Ly49i1",
                    "NOD.Ly49m",
                    "NOD.Ly49u",
                    "129.Ly49g",
                    "129.Ly49o",
                    "129.Ly49p",
                    "129.Ly49v",
                    "B6.Ly49a",
                    "B6.Ly49d",
                    "B6.Ly49f",
                    "B6.Ly49g",
                    "NOD.Ly49a",
                    "NOD.Ly49e",
                    "NOD.Ly49g2",
                    "NOD.Ly49p1",
                    "NOD.Ly49p3")

e4map1.Ly49c <- readRDS("data/e4map1_B6.Ly49c.Rds")
e4map1.Ly49c$strand <- "+"

abao.cage <- aba.subset("figure_map18/abaos/full_v2/abao.Rds", 
                        smart.cut.df = e4map1.Ly49c,
                        regionIDs = paste0(Ly49.cage.hq.3, "$"),
                        new.work.dir = "figure_map18/abaos/e4map1_cage/",
                        threads = 12)

# manual inspection: intron alignment is not perfect, but should be good enough.
# manual inspection: exons and MAPs are well aligned.

# add and write annotations:
abao.cage <- aba.add.track(abao.cage, 
                           track.file = "data/anno_B6.Ly49c_Ly49h_cage.bed",
                           track.regionID.regex = "B6", track.name = "bed",
                           prepend.regionID = T, threads = 1, track.type = "bed")

aba.write.track(abao = abao.cage, consensus.name = "cons_N_0.5",
                write.consensus = T, consensus.add.N = 10000,
                fill.gap = T, broadcast = T, 
                track.name = "bed", track.type = "bed", 
                out.dir = "figure_map18/pileup/cage_v1/bed/")

# add and write nanoCAGE signal tracks:
#### BEGIN: this part was run again with "data/cage_track_df.tsv" for the version with MAPQ==255 filter
abao.cage <- aba.add.track(abao = abao.cage, track.df = "data/noMapq_track_df.tsv",
                           bTrack.df.regionID.regex = T,
                           prepend.regionID = T, threads = 1, track.type = "bdg")

cluster.smart.cut <- utilsFanc::import.bed.fanc("figure_map18/pileup/cage_v1/bed/bed.bed") %>% 
  dplyr::rename(int.name = forth) %>%
  dplyr::mutate(chr = "pos.out") %>% 
  dplyr::filter(int.name %in% c("exon2", "exon1", "RMER5", "MAP1")) %>% 
  dplyr::mutate(end = case_when(int.name == "exon1" ~ as.integer(end + 100), TRUE ~ end))
cluster.smart.cut


aba.write.track(abao = abao.cage,  consensus.name = "cons_N_0.5", track.name = "bdg",
                track.type = "bdg", 
                track.order = function(...) aba.order.by.cat.cage(cluster.blocks = cluster.smart.cut[1:3, ], ...),
                write.consensus = F, normalize = F,
                out.dir = "figure_map18/pileup/cage_v1/cage/", push.2.igv = T)
# note: in IGV view, the "minimum" windowing function was chosen because IGV by default makes noise stand out

#### END.

# plot zoomed-in base pair level heatmaps:
t.f.plot.hm <- function(abao, smart.cut.df, plot.out,
                        width = 10, height = 8, use.max = F, add.nuc = F, 
                        noMapq) {
  bed.file <- "data/anno_B6.Ly49c_Ly49h_cage.bed"
  bdg.df <- paste0("data/", ifelse(noMapq, "noMapq", "cage"), "_track_df.tsv")
  abao <- aba.add.track(abao = abao, track.df = bdg.df,
                        bTrack.df.regionID.regex = T,
                        prepend.regionID = T, threads = 1, track.type = "bdg")
  
  
  track.order <- aba.write.track(abao = abao,  consensus.name = "cons_N_0.5", track.name = "bdg",
                                 track.type = "bdg", 
                                 track.order = function(...) aba.order.by.cat.cage(cluster.blocks = cluster.smart.cut[1:3, ], ...),
                                 normalize = F, order.only = T, out.dir = tempdir())
  
  p1 <- aba.plot.hm(abao = abao, track.type = "bdg", abs = T, 
                    scale.row = F, normalize.row = T, normalize.to.max = use.max,
                    smart.cut.df = smart.cut.df,
                    tracks.include = rev(track.order), 
                    use.order = rev(track.order),
                    cluster.tracks = F, remove.zero.tracks = T,
                    height = 8, use.mafft.order = F,
                    add.nucleotide = add.nuc, project.x.to = "B6.Ly49h")
  
  abao <- aba.add.track(abao, 
                        track.file = bed.file,
                        track.regionID.regex = "B6", track.name = "bed",
                        prepend.regionID = T, threads = 1, track.type = "bed")
  
  p2 <- aba.plot.hm(abao = abao, track.type = "bed",
                    smart.cut.df = smart.cut.df,
                    add.nucleotide = F, broadcast = T, fill.gap = T, bed.same.color = T,
                    project.x.to = "B6.Ly49h"
  )
  p2 <- p2 + theme(legend.position = "bottom")
  
  p3 <- cowplot::plot_grid(p1, p2, align = "v", ncol = 1, rel_heights = c(10,2))  
  dir.create(dirname(plot.out), recursive = T, showWarnings = F)
  ggsave(plot.out, width = width, height = height, dpi = 100)
  return()
}

hm.smart.cut <- read.table("data/anno_B6.Ly49c_Ly49h_cage.bed") %>% 
  `colnames<-`(c("chr", "start", "end", "int.name")) %>% 
  filter(int.name %in% c("exon2", "exon1", "RMER5", "MAP1")) %>% 
  mutate(regionID = case_when(int.name %in% c("exon2", "exon1") ~ "B6.Ly49h",
                              int.name == "MAP1" ~ "B6.Ly49c",
                              int.name == "RMER5" ~ "B6.Ly49i")) %>% 
  mutate(buffer.left = 0, buffer.right = 0)


hm.wMAP1 <- hm.smart.cut %>% 
  mutate(buffer.left = case_when(int.name == "exon2" ~ -20,
                                 int.name == "exon1" ~ -40,
                                 int.name == "RMER5" ~ -80,
                                 int.name == "MAP1" ~ -30)) %>% 
  mutate(buffer.right = case_when(int.name == "exon2" ~ 75,
                                  int.name == "exon1" ~ 90,
                                  int.name == "RMER5" ~ -320,
                                  int.name == "MAP1" ~ -150))

t.f.plot.hm(abao = abao.cage, smart.cut.df = hm.wMAP1, 
            plot.out = "figure_map18/pileup/cage_v1/hm/hm_w_MAP1_max_mapq.pdf",
            width = 20, use.max = T, add.nuc = F, noMapq = F)

# use noMapq = T to generate the version without noMAPQ filter 
