suppressMessages(source("~/R_packages/common/base.R"))
suppressMessages(library(CAGEr))
suppressMessages(library(Rsamtools))
suppressMessages(library(SummarizedExperiment))

for (pkg in c("cageFanc", "bamFanc", "abaFanc", "liteRnaSeqFanc")) {
  # my packages that can be obtained from github
  suppressMessages({
    R.files <- Sys.glob(paste0("~/R_packages/", pkg,"/R/*.R"))
    # modify the locations according to your computer.
    for (R.file in R.files) {
      source(R.file)
    }
  })
}

nkc.shift <- 128534000

strict.filter <- function(in.bams, general.cutoff, threads.master = 3, threads.sub = 4,
                          noMismatch = T, scfilter = T) {
  mclapply(in.bams, function(bam) {
    bam.noMismatch <- sub(".bam$", "_noMismatch.bam", bam)
    bam.scfilter <- sub(".bam$", "_scfilter.bam", bam.noMismatch)
    if (noMismatch == T) {
      stream.bam.core(bam.file = bam, fields = ALL.FIELDS, tags = CAGE.TAGS,header.file = NULL,
                      chunk.size = 30000, chunk.call.back = bam.filter.chunk,
                      chunk.call.back.param.list = list(outfile = bam.noMismatch,
                                                        threads = threads.sub,
                                                        filter.fun = bam.filter.no.mismatch,
                                                        filter.fun.params = list(tlen.cutoff = 10000)))
    }
    if (scfilter == T) {
      stream.bam.core(bam.file = bam.noMismatch, fields = ALL.FIELDS, tags = CAGE.TAGS,
                      chunk.size = 30000, chunk.call.back = bam.filter.chunk,
                      chunk.call.back.param.list = list(outfile = bam.scfilter,
                                                        # junkfile = "test/sth/mini_junk.bam",
                                                        # remove.mate = T,
                                                        filter.fun = bam.filter.softclip,
                                                        filter.fun.params = list(general.cutoff = general.cutoff,
                                                                                 threads = threads.sub)),
                      header.file = NULL)
    }
    
  }, mc.cores = threads.master, mc.cleanup = T)
}



se.get.info <- function(se, assay, annotate.se = F, 
                        # to annotate RangedSummarizedExperiment
                        genome,
                        rowData.cols = NULL, colData.cols = NULL, # use all when NULL
                        rowData.filter.col = NULL, colData.filter.col = NULL, # use rownames and colnames when NULL
                        rows.include = NULL, cols.include = NULL, # include everything when NULL
                        use.regex = T, use.Ly49 = T, value.name = "value", 
                        return.se = F) {
  if (annotate.se == T) {
    se <- utilsFanc::gr.fast.annot(obj = se, genome = genome, anno.cols = rowData.cols)
  }
  
  if (!is.null(rowData.cols)) {
    rowData(se) <- rowData(se) %>% .[, colnames(.) %in% c(rowData.cols, rowData.filter.col), drop = F]
  }
  if (!is.null(colData.cols)) {
    colData(se) <- colData(se) %>% .[, colnames(.) %in% c(colData.cols, colData.filter.col), drop = F]
  }
  
  if (is.null(rowData.filter.col)) {
    rowData.filter.col <- "feature"
    if (is.null(rownames(se))) {
      rownames(se) <- 1:nrow(se)
    }
    rowData(se)$feature <- rownames(se)
  }
  if (is.null(colData.filter.col)) {
    colData.filter.col <- "sample"
    colData(se)$sample <- colnames(se)
  }
  if (! rowData.filter.col %in% colnames(rowData(se))) {
    stop("! rowData.filter.col %in% colnames(rowData(se))")
  }
  
  if (! colData.filter.col %in% colnames(colData(se))) {
    stop("! colData.filter.col %in% colnames(colData(se))")
  }
  
  if (!is.null(rows.include)) {
    if (use.regex == T) {
      include.regex <- paste0(rows.include, collapse = "|")
      rows.include <- rowData(se)[, rowData.filter.col] %>% .[grepl(include.regex,.)]
    }
    se <- se %>% .[rowData(.)[, rowData.filter.col] %in% rows.include, ]
    if (nrow(se) < 1) {
      stop("nrow(se) < 1")
    }
  }
  
  if (!is.null(cols.include)) {
    if (use.regex == T) {
      include.regex <- paste0(cols.include, collapse = "|")
      cols.include <- colData(se)[, colData.filter.col] %>% .[grepl(include.regex,.)]
    }
    se <- se %>% .[, colData(.)[, colData.filter.col] %in% cols.include]
    if (ncol(se) < 1) {
      stop("ncol(se) < 1")
    }
  }
  if (use.Ly49 == T) {
    rowData(se)[, rowData.filter.col] <- rowData(se)[, rowData.filter.col] %>%
      utilsFanc::klra.ly49.translate()
  }
  # se ready
  if (return.se == T) {
    if (is.null(rowData.filter.col)) {
      rowData[, rowData.filter.col] <- NULL
    }
    if (is.null(colData.filter.col)) {
      colData[, colData.filter.col] <- NULL
    }
    
    return(se)
  }
  if ("RangedSummarizedExperiment" %in% class(se)) {
    rowData(se)$locus <- utilsFanc::gr.get.loci(gr = rowRanges(se), keep.strand = T)
  }
  df <- cbind(rowData(se), assay(se, assay)) %>% as.data.frame()
  df <- df %>% reshape2::melt(id.vars = colnames(rowData(se)), 
                              measure.vars = colnames(se),
                              value.name = value.name, 
                              variable.name = "sample")
  df <- left_join(df, colData(se) %>% as.data.frame() %>% mutate(., sample = rownames(.)))
  return(df)
}

