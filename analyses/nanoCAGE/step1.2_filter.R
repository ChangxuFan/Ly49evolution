source("objects.R")

### MAPQ255 (default from cageFast.R)
bams <- Sys.glob("sth/star_aligned/B6*/*.out.umi.dedup.filtered.bam") %>%
  c(Sys.glob("sth/*/star_aligned/*/*.out.umi.dedup.filtered.bam"))

# remove reads with mismatches and extensive softclipping
strict.filter(in.bams = bams, general.cutoff = 3, threads.master = 1, threads.sub = 10)
# cap filter
B6.bams <- Sys.glob("~/4dn/nk/cage/nk/sth/star_aligned/*/*noMismatch_scfilter.bam")
nod129.bams <- c("~/4dn/nk/cage/nk/sth/129/star_aligned/m129_mRNA/m129_mRNA_Aligned.sortedByCoord.out.umi.dedup.filtered_noMismatch_scfilter.bam",
                 "~/4dn/nk/cage/nk/sth/NOD/star_aligned/NOD_mRNA/NOD_mRNA_Aligned.sortedByCoord.out.umi.dedup.filtered_noMismatch_scfilter.bam")

cap.filter.bam(B6.bams, npar = 4)
cap.filter.bam(nod129.bams, npar = 2)

### END: MAPQ255 (default from cageFast.R)

### no MAPQ filter
bams <- c(Sys.glob("sth/star_aligned/B6*/B6*dedup.bam"),
          "sth/129/star_aligned/m129_mRNA/m129_mRNA_Aligned.sortedByCoord.out.umi.dedup.bam", 
          "sth/NOD/star_aligned/NOD_mRNA/NOD_mRNA_Aligned.sortedByCoord.out.umi.dedup.bam")
out.dir <- "noMapq_filter/sth/bams/"

utilsFanc::safelapply(bams, function(bam) {
  sample <- basename(dirname(bam))
  out.bam <- paste0(out.dir, "/", sample, "/", sample, ".umi.dedup.noMapq.bam")
  samtools.filter.juheon(in.bam = bam, out.bam = out.bam, 
                         no.mapq = T,
                         thread = 3, run = T)
  return(out.bam)
}, threads = 6)

bams <- Sys.glob("noMapq_filter/sth/bams/*/*dedup.noMapq.bam")
strict.filter(in.bams = bams, general.cutoff = 3, threads.master = 6, threads.sub = 3)

bams <- Sys.glob("noMapq_filter/sth/bams/*/*noMismatch_scfilter.bam")
cap.filter.bam(bams, npar = 6)

filtered.bam <- bams %>% 
  sub(".bam$", "_capFilter.bam", .)

if(any(!file.exists(filtered.bam))) {
  stop()
}

### END: no MAPQ filter