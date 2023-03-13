source("~/R_packages/common/base.R")
source("~/R_packages/cageFanc/R/values.R")
source("~/R_packages/cageFanc/R/star_wrapper.R")

fastqs <- Sys.glob("fastq/*fastq.gz")
# star genome index was made using:
# STAR --runThreadN 12 --runMode genomeGenerate \
# --genomeDir ~/genomes/mm10/STAR_gencode_vM24 --genomeFastaFiles ~/genomes/mm10/mm10.fa \
# --sjdbGTFfile ~/genomes/mm10/gencode/gencode.vM24.annotation.gtf --sjdbOverhang 100
star.fanc.2(fastqs = fastqs, 
            genome.index = "~/genomes/mm10/STAR_gencode_vM24/",
            threads.sample = 4, threads.each = 4,
            outdir = "star_aligned/")
