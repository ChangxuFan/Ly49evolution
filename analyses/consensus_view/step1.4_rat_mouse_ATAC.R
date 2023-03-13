source("objects.R")
aba.all <- readRDS("aba_df_rat_mouse.Rds")
aba.create(aba.df = aba.all, df.is.bed = F, work.dir = "pileup/abaos/e1m1_add_v0/",
           mask.lower.regex = c("B6.Ly49q", "B6.Ly49alpha"), threads = 12)

# examine MAP8 alignment
map8.df <- data.frame(chr = "chr6", start = 130341136, end = 130341386,
                      regionID = "B6.Ly49c", genome = "mm10")

trash <- aba.subset("pileup/abaos/e1m1_add_v0/abao.Rds", smart.cut.df = map8.df,
                    new.work.dir = "pileup/abaos/e1m1_add_v0_map8/")
# visually alignment looks good

# trim off the unaligned fractions at the 2 ends
trash <- aba.subset("pileup/abaos/e1m1_add_v0/abao.Rds", 
                    smart.cut.df = aba.all[aba.all$regionID == "B6.Ly49c",],
                    new.work.dir = "pileup/abaos/e1m1_add/",
                    mask.lower.regex = c("B6.Ly49q", "B6.Ly49alpha"), threads = 12)


######### now we plot tracks~
abao.e1m1 <- readRDS("pileup/abaos/e1m1_add/abao.Rds")
abao.e1m1 <- aba.add.consensus(abao.e1m1)

abao.e1m1 <- aba.add.track(abao.e1m1, 
                           track.file = "data/anno_B6.Ly49c_simple.bed",
                           track.name = "bed",
                           track.regionID.regex = "B6",
                           prepend.regionID = T, threads = 1, track.type = "bed")

aba.write.track(abao = abao.e1m1, consensus.name = "cons_N_0.5",
                write.consensus = T, consensus.add.N = 10000,
                fill.gap = T, broadcast = T, 
                track.name = "bed", track.type = "bed", 
                out.dir = "pileup/v4_add/bed/")


track.dfs <- lapply(list(q0 = 0, q8 = 8), function(mapq) {
  lapply(list(NK = "NK"), function(type) {
    if (mapq == 0) {
      mapq <- ""
    } else {
      mapq <- "_mapq8"
    }
    df <- data.frame(regionID = c("B6", "rn7"),
                     track.name = "atac", track.type = "bdg",
                     track.file = c(
                       paste0("/bar/cfan/4dn/nk/fanc/atac2/sth/snakeMAPQ1/methylQA/Dpos_rep1_AS0", mapq, ".bigWig"),
                       paste0("~/4dn/nk/rat/mcw/sth/snakeBothRep/methylQA/rat", type, "_rep1_AS0", mapq, ".bigWig")
                     ))
    return(df)
  })
})

abao.e1m1 <- aba.add.track(abao = abao.e1m1, track.df = track.dfs$q0$NK,
                           bTrack.df.regionID.regex = T,
                           prepend.regionID = T, threads = 1, track.type = "bdg")

aba.write.track(abao = abao.e1m1, consensus.name = "cons_N_0.5", track.name = "bdg",
                track.type = "bdg", track.order = "aln", write.consensus = F, 
                normalize = T, add.constant = 1.5,
                out.dir = "pileup/v4_add/", push.2.igv = T)

# to get the view for the MAPQ >= 8 filtered version: change track.df to track.dfs$q8$NK
# manually remove from alignment pileup on IGV:
# s7/s8: assembly issues (detailed in Methods)
#########