source("objects.R")
Klra3.gtf <- rtracklayer::import("../../annotations/gene_annotations/Ly49_mm10.gtf") %>% 
  .[.$gene_name == "B6.Ly49c"]
ma.gtf <- rtracklayer::import("../../annotations/gene_annotations/mesAur_2.0_genes.gtf")

# first align mesAur to Ly49c (Klra3) to check alignment quality
t <- lapply(1:2, function(intron.id) {
  gtfs <- list(ma = ma.gtf, Klra3 = Klra3.gtf)
  aba.df <- lapply(names(gtfs), function(name) {
    gtf <- gtfs[[name]]
    df <- gtf %>% utilsFanc::gr2df() %>% 
      dplyr::filter(exon_number == intron.id, type == "intron") %>% 
      dplyr::select(chr, start, end, strand, gene_name) %>% 
      dplyr::rename(regionID = gene_name)
    if (name == "ma") {
      df$genome <- "mesAur2.0"
      if (intron.id == 2) {
        df <- df[df$regionID != "ma.Ly49ma5",]
      }
    } else if (name == "Klra3") {
      df$genome <- "mm10"
    } else {
      stop()
    }
    df$strand <- ifelse(df$strand == "+", "-", "+")
    return(df)
  }) %>% do.call(rbind, .)
  trash <- aba.create(aba.df = aba.df,
                      work.dir = paste0("mesAur/abaos_intron/intron", intron.id, "_ma_Klra3/"),
                      df.is.bed = F)
  return(aba.df)
})

# END: first align mesAur to Ly49c (Klra3) to check alignment quality

# select regions used for gene tree construction (refer to Methods)
# make sure mesAur sequences can align to Ly49c
Klra3.cuts <- lapply(1:2, function(intron.id) {
  df <- readRDS(paste0("mesAur/abaos_intron/intron", intron.id, "_ma_Klra3/abao.Rds"))@ori.df %>% 
    dplyr::filter(regionID == "B6.Ly49c")
  if (intron.id == 1) {
    df$start <- df$start + 200
    df$end <- df$end - 200
  } else if (intron.id == 2) {
    df$end <- df$start + 913
    df$start <- df$start + 200
  }
  return(df)
})

names(Klra3.cuts) <- paste0("intron", 1:2)

lapply(1:2, function(intron.id) {
  trash <- aba.subset(paste0("mesAur/abaos_intron/intron", intron.id, "_ma_Klra3/abao.Rds"),
                      smart.cut.df = Klra3.cuts[[intron.id]],
                      new.work.dir = paste0("mesAur/abaos_intron/intron", intron.id, "_ma_Klra3_cut/"))
  return()
})

# END: select regions used for gene tree construction (refer to Methods)

# now we align mesAur Ly49 sequences to all mouse/rat Ly49 sequences.
mr.dfs <- lapply(1:2, function(intron.id) {
  df <- aba.subset(paste0("mesAur/abaos_intron/mr_i2i1/i", intron.id, "/abao.Rds"),
                   smart.cut.df = Klra3.cuts[[intron.id]], ori.df.only = T)
  return(df)
})

ma.use <- paste0("ma.Ly49ma", c(1,3,6,7))

lapply(1:2, function(intron.id) {
  ma.df <- readRDS(paste0("mesAur/abaos_intron/intron", intron.id, "_ma_Klra3_cut/abao.Rds"))@ori.df %>% 
    dplyr::filter(genome == "mesAur2.0", regionID %in% ma.use)
  aba.df <- rbind(ma.df, mr.dfs[[intron.id]])
  trash <- aba.create(aba.df, work.dir = paste0("mesAur/abaos_intron/intron", intron.id, "_mrma/"), df.is.bed = F,
                      threads = 12)
  return()
})

# END: now we align mesAur Ly49 sequences to all mouse/rat Ly49 sequences.
fas <- lapply(1:2, function(intron.id) {
  mb.dir <- paste0("mesAur/abaos_intron/intron", intron.id, "_mrma/cons/cons_N_0.5/mb/")
  fa <- Sys.glob(paste0("mesAur/abaos_intron/intron", intron.id, "_mrma/cons/cons_N_0.5/*_shrink.fa"))
  if (length(fa) < 1) {
    stop()
  }
  return(fa)
}) %>% unlist()
names(fas) <- c("i1", "i2")

# prepare input files for mrBayes
fa.2.nex(fa.vec = fas, out.file = "mesAur/abaos_intron/intron12_mrma_mrbayes/intron12_mrma_shrink.nex", 
         use.partition = T, add.mb.GTR = T)
# mrBayes was run according to Methods
# END: prepare input files for mrBayes

# remove branches with low posterior probability
mb.tree.shorten.branch(tree.in = "mesAur/abaos_intron/intron12_mrma_mrbayes/intron12_mrma_shrink.nex.con.tre")
# END: remove branches with low posterior probability


