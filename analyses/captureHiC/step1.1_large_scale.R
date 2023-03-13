source("objects.R")

Klra4.25kb <- "chr6:130051712-130072863" %>% utilsFanc::loci.2.gr()
nkc.blocks <- c("chr6:129582583-129744372", "chr6:129764477-130015298", "chr6:130105696-130391935") %>% 
  utilsFanc::loci.2.gr()
nkc.blocks$block.name <- c("Nkg2", "Ly49_inactive", "Ly49_active")

hico <- read.hic(.hic.vec = INTER.hic, region1 = REGION.large.v1, observed.or.oe = "observed", 
                 normalization = c("NONE", "VC_SQRT", "VC"), use.1.6 = T,
                 threads.norm = NULL, threads.res = 1, threads.sample = NULL, 
                 colData = DataFrame(Ly49D = c("neg", "pos", "neg", "pos"),
                                     rep = c("rep1", "rep1", "rep2", "rep2")))
dir.create("publication")
saveRDS(hico, "publication/hico_v1.Rds")
# a copy of hico is provided 

ll.m <- list()
ll.m$v6.equal.v1 <- contact.pipe(
  anchor1 = Klra4.25kb, anchor2 = nkc.blocks, resolution = 25000,
  gi.name.col = "anchor2.block.name", 
  intSet = hico$res_25000, assay = "VC", sum.fun = sum,
  out.dir = "publication/large_scale/v6.equal.v1/", root.name = "v6.equal.v1", 
  write.anchors = T)

# make the plot
df <- readRDS("publication/large_scale/v6.equal.v1/v6.equal.v1.Rds")$contact.stat$stat.melt

df <- df %>% dplyr::filter(int.name != "Nkg2") %>% 
  dplyr::mutate(D = ifelse(grepl("neg", sample), "Ly49D-", "Ly49D+"),
                rep = paste0("rep", stringr::str_extract(sample, "\\d")),
                block = ifelse(int.name == "Ly49_inactive", "ncNK_Ly49", "cNK_Ly49")) %>% 
  dplyr::group_by(D, rep) %>% 
  dplyr::mutate(rel.freq = int.freq/sum(int.freq)) %>% 
  dplyr::ungroup() %>% as.data.frame() %>% 
  dplyr::arrange(D, rep, block) %>% 
  dplyr::filter(block == "cNK_Ly49") %>% 
  dplyr::select(D, rep, block, rel.freq)

df$D <- df$D %>% sub("Ly49", "", .)
df$D.color <- df$D
df$D <- factor(df$D, levels = c("D-", "D+"))
df$D.color <- factor(df$D.color, levels = c("D+", "D-"))
p <- utilsFanc::barplot.pub.3(
  df = df, x = "D", color.by = "D.color", y = "rel.freq",
  palette.fc = "red_green", bar.width = 0.8, spread.width = 0.8
) %>% 
  utilsFanc::theme.fc.1(italic.x = F) +
  ggsave("publication/large_scale/large_scale.pdf", device = cairo_pdf, 
         width = 0.75, height = 1, dpi = 300)
