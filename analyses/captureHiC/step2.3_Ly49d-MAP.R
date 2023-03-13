rm(list = ls())
source("objects.R")
enhancer.summits <- readRDS("publication/enhancer.summits.Rds")
enhancer.bins <- enhancer.summits %>% utilsFanc::gr.fit2bin(bin.size = 5000, expand = T)
enhancer.bins <- gr.use.ly49(enhancer.bins)
Klra4.enhancer.gi <- GInteractions(anchor1 = rep(1, length(enhancer.bins)), 
                                   anchor2 = 2:(length(enhancer.bins) + 1), 
                                   regions = c(KLRA4.bins$res_5000$promoter, enhancer.bins))
regions(Klra4.enhancer.gi)
Klra4.enhancer.gi$int.name <- paste0(anchors(Klra4.enhancer.gi, "first")$forth, "_", 
                                     anchors(Klra4.enhancer.gi, "second")$forth)
hico.list <- readRDS("publication/hico.list.v1.Rds")
names(hico.list) <- sub("30", "_30", names(hico.list))
normo <- readRDS("publication/SIPG/coverage/normo_v1.Rds")
trash <- SIPG.pipe(fg.gi = Klra4.enhancer.gi, hico.list = hico.list, 
                   normo.list = normo, fg.ext = 15000, 
                   fg.id.col = "int.name", distance.wobble = 5000, resolutions = 5000, 
                   mapqs = c("inter.hic", "inter_30.hic"), norms = c("VC"),
                   out.dir = "publication/SIPG/Klra4_enhancer_grid_ext_15000/", 
                   color.range = NULL, center.piece = c(4, 4), debug = T,
                   plot.only = F, plot.mat.grid = F,
                   plot.bar = F, replot.collapse = F, max.quantile = 0.999,
                   sample.coldata = SAMPLE.COLDATA.hic
)

#### generate figures for step2.2 and 2.3
extrafont::loadfonts()
## to the reader: you would need to install the Arial font

lapply(c("map", "Ly49d"), function(x) {
  if (x == "map") {
    data <- readRDS("publication/SIPG/EL4_grid_ext_15000/inter.hic..res_5000..VC/result.Rds")$center.df
  } else {
    data <- readRDS("publication/SIPG/Klra4_enhancer_grid_ext_15000/inter.hic..res_5000..VC/result.Rds")$center.df
  }
  
  data$sample <- data$sample %>% stringr::str_extract("NK\\d.D.") # %>% 
  # sub("NK", "", .) %>% sub(".D", "",.) %>% sub("(\\d)([np])", "\\2\n\\1")
  data$sample[is.na(data$sample)] <- "EL4"
  levels <- c("NK1-Dn", "NK1-Dp", "NK2-Dn", "NK2-Dp")
  if (x == "map") {
    levels <- c(levels, "EL4")
  }
  data$sample <- factor(
    data$sample, levels = levels)
  data$type <- factor(data$type, c("fg", "bg"))
  plot.dir <- "publication/SIPG/plots/"
  dir.create(plot.dir, showWarnings = F, recursive = T)
  plot.out <- paste0(plot.dir, "/", x, "_map.pdf")
  p <- utilsFanc::barplot.pub.3(df = data, x = "sample", y = "center.pct",
                     color.by = "type", jitter.dots = T,
                     palette.fc = "red_green",
                     bar.width = 0.8, dodge.width = 0.9, bar.line.size = 0.5, 
                     error.bar.width = 0.4, error.bar.line.size = 0.3,
                     add.pval = T, pval.group.1 = "bg", pval.group.2 = "fg", 
                     pval.same.y = T,
                     pt.size = 0.5) %>% 
    utilsFanc::theme.fc.1(italic.x = F) +
    ggsave(plot.out, device = "pdf", 
           width = 1.75, height = 1.0, dpi = 300)
  embedFonts(plot.out)
})

# map-map barplot (NK vs EL4) and Ly49d promoter-map barplot (Ly49D+ vs Ly49D-)
