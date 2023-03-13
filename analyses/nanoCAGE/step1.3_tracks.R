source("objects.R")

### MAPQ255 (default from cageFast.R)
B6.bams <- Sys.glob("~/4dn/nk/cage/nk/sth/star_aligned/*/*noMismatch_scfilter.bam")
nod129.bams <- c("~/4dn/nk/cage/nk/sth/129/star_aligned/m129_mRNA/m129_mRNA_Aligned.sortedByCoord.out.umi.dedup.filtered_noMismatch_scfilter.bam",
                 "~/4dn/nk/cage/nk/sth/NOD/star_aligned/NOD_mRNA/NOD_mRNA_Aligned.sortedByCoord.out.umi.dedup.filtered_noMismatch_scfilter.bam")

filtered.bam <- c(B6.bams, nod129.bams) %>% 
  sub(".bam$", "_capFilter.bam", .)
file.exists(filtered.bam)

filtered.bam.df <- data.frame(bam = filtered.bam, 
                              sample = sub("_Aligned.+$", "", basename(filtered.bam)),
                              genome = c(rep("mm10", 4), NA, NA),
                              genome.name = c(rep(NA, 4), "BSgenome.Mmusculus.fanc.129S1Ly49S6",
                                              "BSgenome.Mmusculus.fanc.NODLy49fix"))

ce.l.filter <- bam2tss.fanc(bam.df = filtered.bam.df, out.dir = "filter/sth/tracks_capFilter/", threads.master = 1,
                            threads.sub = 8)
### END: MAPQ255 (default from cageFast.R)

### no MAPQ filter
filtered.bam.df <- data.frame(bam = filtered.bam, 
                              sample = sub(".umi.+$", "", basename(filtered.bam)),
                              genome = c(rep("mm10", 4), NA, NA),
                              genome.name = c(rep(NA, 4), "BSgenome.Mmusculus.fanc.129S1Ly49S6",
                                              "BSgenome.Mmusculus.fanc.NODLy49fix"))

ce.l.filter <- bam2tss.fanc(bam.df = filtered.bam.df, out.dir = "filter/sth/tracks_noMapq_capFilter/", threads.master = 1,
                            threads.sub = 8)

### END: no MAPQ filter