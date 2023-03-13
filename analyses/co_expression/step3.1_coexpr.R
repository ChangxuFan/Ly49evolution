rm(list = ls())
source("objects.R")
soi <- readRDS("vivier_2022/sth/merge_v3_MM/soi.Rds")
soi <- soi[, soi$celltype == "NK"]
# this file is provided for your convenience.
Klra.expressed <- paste0("Klra", c(6, 4, 8, "14-ps", 9, 7, 10, "13-ps", 3, 1))
# please note that this function has a lot of code initially written 
# for other purposes... 
# method = "pos.ratio" instructs the code to perform the odds ratio analyses 
# used in this paper.
coexpr.grid(so = soi, assay = "RNA", slot = "counts", 
            genes = Klra.expressed,  
            gene.alias = "translate_strip", 
            use.gene.order = T,
            method = "pos.ratio", pos.ratio.log2 = F,
            hm.values = c(0.5, 1, 3), gene.name.fontsize = 6,
            out.dir = "vivier_2022/coexpr/pos.ratio/", pub = T)
