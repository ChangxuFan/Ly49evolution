source("objects.R")

hico.list <- lapply(list(INTER.hic, INTER30.hic), function(.hic.vec) {
  hico <- read.hic(.hic.vec = .hic.vec, region1 = REGION.large.v1, observed.or.oe = "observed", 
                   normalization = c("NONE", "VC_SQRT", "VC"), use.1.6 = T,
                   threads.norm = NULL, threads.res = 1, threads.sample = NULL, 
                   colData = DataFrame(Ly49D = c("neg", "pos", "neg", "pos"),
                                       rep = c("rep1", "rep1", "rep2", "rep2")))
  return(hico)
})

names(hico.list) <- c("inter.hic", "inter30.hic")
saveRDS(hico.list, "publication/hico.list.v1.Rds")
# a copy is provided 

Klra4.body <- "chr6:130044017-130067271" %>% utilsFanc::loci.2.gr()
partner1.body <- "chr6:130291161-130306432" %>% utilsFanc::loci.2.gr()
width(Klra4.body)
Klra4.exp20k <- Klra4.body + 20000
partner1.exp20k <- partner1.body + 20000

anchor.list <- list(v1 = list())
anchor.list$v1$a1 <- c(Klra4.body, Klra4.exp20k)
anchor.list$v1$a2 <- c(partner1.body, partner1.exp20k)
anchor.list$v1$gi <- GInteractions(anchor.list$v1$a1, anchor.list$v1$a2)
anchor.list$v1$gi$int.name <- c("Ly49d_partner1", "peripheral")
anchor.list$v1$gi
dot.m <- list()
dot.m$p1 <- contact.pipe.m(
  hico.by.mapq = hico.list, norms = c("VC", "NONE", "VC_SQRT"), 
  gi = anchor.list$v1$gi, gi.name.col = "int.name", sum.fun = mean, 
  resolution = 5000, out.dir = "publication/dot/v2_grid_5K_repeat_mean/", 
  do.donut = T, donut.center = "Ly49d_partner1", donut.peri = "peripheral",
  coldata = SAMPLE.COLDATA)

#### now we move on to the other partner:

Klra4.body <- "chr6:130044017-130067271" %>% utilsFanc::loci.2.gr()
partner2.body <- "chr6:130149106-130160748" %>% utilsFanc::loci.2.gr()

al.2 <- list()
Klra4.body <- "chr6:130044017-130067271" %>% utilsFanc::loci.2.gr()
partner2.body <- "chr6:130149106-130160748" %>% utilsFanc::loci.2.gr()

al.2$v1$a1 <- c(Klra4.body, Klra4.body + 20000)
al.2$v1$a2 <- c(partner2.body, partner2.body + 20000)

al.2$v1$gi <- GInteractions(al.2$v1$a1, al.2$v1$a2)
al.2$v1$gi$int.name <- c("Ly49d_partner2", "peripheral")
al.2$v1$gi

dot.m$p2 <- contact.pipe.m(
  hico.by.mapq = hico.list, norms = c("NONE", "VC_SQRT", "VC"), 
  gi = al.2$v1$gi, gi.name.col = "int.name", sum.fun = mean, 
  resolution = 5000, out.dir = "publication/dot2/v1.grid.5k_mean/", 
  do.donut = T, donut.center = "Ly49d_partner2", donut.peri = "peripheral",
  coldata = SAMPLE.COLDATA)

dir.create("publication/dots/")
saveRDS(dot.m, "publication/dots/dot.m.Rds")
#### now generate figures:

lapply(c("partner1", "partner2"), function(x) {
  if (x == "partner1") {
    data <- readRDS("publication/dot/v2_grid_5K_repeat_mean/grid.Rds")
  } else {
    data <- readRDS("publication/dot2/v1.grid.5k_mean/grid.Rds")
  }
  loop.name <- paste0("Ly49d_", x)
  data <- data$inter.hic..VC$contact.stat$stat.melt
  data$Ly49D[data$Ly49D == "neg"] <- "D-"
  data$Ly49D[data$Ly49D == "pos"] <- "D+"
  data$Ly49D <- factor(data$Ly49D, levels = c("D-", "D+"))
  data$int.name <- factor(data$int.name, levels = c(loop.name, "peripheral..donut"))
  print(data)
  p <- utilsFanc::barplot.pub.3(
    df = data, x = "Ly49D", y = "int.freq",
    color.by = "int.name", palette.fc = "red_green",
    bar.width = 0.6, dodge.width = 0.7, bar.line.size = 0.6, 
    spread.bin.size = 0.001,
    add.pval = T, pval.use.star = T, 
    pval.group.1 = "peripheral..donut", 
    pval.group.2 = loop.name,
    pval.bar.y.nudge = 0.12, 
    pval.text.y.nudge = 0.18) %>% 
    utilsFanc::theme.fc.1(italic.x = F) +
    theme(aspect.ratio = 1) +
    ggsave(paste0("publication/dots/d", x, "_loop.pdf"), device = cairo_pdf, 
           width = 1, height = 1, dpi = 300)
  return()
})

####


