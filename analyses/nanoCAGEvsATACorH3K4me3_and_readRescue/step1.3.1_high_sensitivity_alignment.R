source("objects.R")
fastqs <- Sys.glob("/PATH/to/your/cutadapt/*.fastq.gz")

names(fastqs) <- fastqs %>% stringr::str_extract("\\d-D[(neg)(pos)]") %>% tolower() %>% 
  sub("(\\d)-(.+)", "\\2\\1", .)

bam.dir <- "nkc_bowtie2/sth/atac_150_nkc"
dir.create(bam.dir, showWarnings = F, recursive = T)
fastqs %>% split(f = names(fastqs)) %>% mclapply(function(fastq) {
  out.bam <- bam.dir %>% paste0("/", names(fastq)[1], ".bam")
  out.log <- bam.dir %>% paste0("/", names(fastq)[1], ".log")
  liteRnaSeqFanc::bowtie2.wrapper(
    fastq.or.pair = fastq, genome.index = "~/genomes/nkc/bowtie2/nkc",
    out.bam = out.bam, 
    threads = 3, a = T, X = 5000, mm = T, report.unaligned = F, 
    add.params = " --no-mixed ", log.file = out.log, 
    preset = " --very-sensitive ", run = T,
    bowtie2 = "/opt/apps/bowtie2/2.3.4.1/bowtie2")
}, mc.cores = 6, mc.cleanup = T)

# genome index: bowtie2 index for mm10:chr6:128534000-131409000
