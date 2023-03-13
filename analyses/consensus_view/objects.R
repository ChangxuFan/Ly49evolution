source("~/R_packages/common/base.R")
source("~/R_packages/common/bioc.R")
suppressMessages(library(grid))
suppressMessages(library(gridExtra))
suppressMessages(library(ComplexHeatmap))
suppressMessages(library(ape))
suppressMessages(library(phangorn))
suppressMessages(library(ggtree))
suppressMessages(library(treeio))
suppressMessages(library(ggrepel))
suppressMessages(library(ggpubr))
suppressMessages(library(circlize))
suppressMessages(library(abaFanc2))
suppressMessages(library(liteRnaSeqFanc))

PATHs <- Sys.getenv("PATH") %>% strsplit(split = ":") %>% unlist()
if (! BEDTOOLS.DIR %in% PATHs) {
  PATHs <- c(PATHs, BEDTOOLS.DIR)
  Sys.setenv(PATH = PATHs %>% paste0(collapse = ":"))
}

if (! HOMER.BIN %in% PATHs) {
  PATHs <- c(PATHs, HOMER.BIN)
  Sys.setenv(PATH = PATHs %>% paste0(collapse = ":"))
}

grep.acti <- paste0("Ly49", c("h", "n","k", "m", "d", "p", "pd", "l", "r", "lr", "x", "u", "ui")) %>% 
  paste0(collapse = "|")

grep.inhi <- paste0("Ly49", c("a", "o", "g", "v", "t", "c", "e", "ec", "f", "i", "j", "s") )%>% 
  paste0(collapse = "|")
grep.list <- list(acti = grep.acti, inhi = grep.inhi)
grep.list$B6.acti <- paste0("B6.Ly49", c("h", "n","k", "m", "d", "x"))
# gene.names <- readLines("/bar/cfan/4dn/nk/evolution/exon/per_exon_align/CDS_2/gene_names.txt")
# length(gene.names)
# acti <- gene.names %>% .[grepl(grep.acti, .)]
# inhi <- gene.names %>% .[grepl(grep.inhi, .)]
# length(acti)
# length(inhi)
# 
gg_color_hue <- utilsFanc::gg_color_hue

bed.2.cat <- function(beds, out.json = NULL) {
  if (is.null(out.json)) {
    if (length(beds) == 1) {
      out.json <- beds %>% sub(".gz$", "",.) %>% paste0(".cat.json")
    }
  }
  jgs <- lapply(beds, function(bed) {
    bed <- bed %>% sub(".gz$", "", .)
    df <- read.table(bed, header = F, sep = "\t", quote = "")
    cats <- df$V4 %>% unique() %>% sample(., length(.), replace = F)
    colors <- gg_color_hue(length(cats))
    cat.list <- lapply(seq_along(cats), function(i) {
      return(list(name = cats[i], color = colors[i]))
    })
    names(cat.list) <- cats
    cat.options <- list(category = cat.list)
    other.options <- list(height = 20, alwaysDrawLabel= T, maxRows= 100, hiddenPixels= 0)
    jg <- list(type = "categorical", 
               name = basename(bed),
               url = utilsFanc::bash2ftp(paste0(bed, ".gz")),
               options = c(other.options, cat.options))
    return(jg)
  })
  json <- jgs %>% jsonlite::toJSON(auto_unbox = T) %>% jsonlite::prettify()
  write(json, out.json)
  cat(utilsFanc::bash2ftp(out.json))
  cat("\n")
  return()
}

ensemblID.2.gene <- function(ensembl_ID, genome) {
  # cow: btaurus; human: hsapiens
  ensembl <- biomaRt::useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL")
  ensembl <- biomaRt::useDataset(paste0(genome, "_gene_ensembl"), mart = ensembl)
  G_list <- biomaRt::getBM(filters = "ensembl_gene_id", attributes = c("ensembl_gene_id","external_gene_name"),
                           values = ensembl_ID, mart= ensembl)
  ori.df <- data.frame(ensembl_gene_id = ensembl_ID)
  df <- left_join(ori.df, G_list)
  if (nrow(df) != nrow(ori.df)) {
    stop("nrow(df) != nrow(ori.df)")
  }
  return(df$external_gene_name)
}

FAANG.rna <- function(files, genome) {
  df <- lapply(names(files), function(sample) {
    df <- read.table(files[[sample]])
    colnames(df) <- c("ensembl_ID", sample)
    return(df)
  }) %>% Reduce(full_join, .)
  rownames(df) <- df$ensembl_ID
  df$gene <- ensemblID.2.gene(df$ensembl_ID, genome = genome)
  df <- df[!is.na(df$gene) & df$gene != "", ]
  df$ensembl_ID <- NULL
  df <- df %>% group_by(gene) %>% 
    summarise_all(sum) %>% ungroup() %>% as.data.frame()
  
  rownames(df) <- df$gene
  df$gene <- NULL
  mat <- df %>% as.matrix()
  mat <- edgeR::cpm(mat)
  return(mat)
}

featureCounts.2.cpm <- function(file) {
  df <- read.table(file, header = T)
  rownames(df) <- df$gene
  df$gene <- NULL
  mat <- as.matrix(df)
  mat <- edgeR::cpm(mat)
  return(mat)
}

expr.bar <- function(mat, genes, plot.out = NULL,
                     sub.width = 10, sub.height = 4, n.col = 1,
                     n.show = NULL) {
  pl <- lapply(genes, function(gene) {
    df <- data.frame(sample = colnames(mat), expr = mat[gene,]) %>% 
      arrange(desc(expr))
    if (is.null(n.show)) {
      n.show <- nrow(df)
    }
    df <- df[1:n.show, ]
    df <- df %>% arrange(expr) %>% 
      dplyr::mutate(sample = factor(sample, levels = sample))
    p <- ggplot(df, aes(x = sample, y = expr)) +
      geom_bar(stat = "identity", alpha = 0.5, fill = "purple4") +
      scale_x_discrete(guide = guide_axis(angle = 90))
    p <- p %>% utilsFanc::theme.fc.1(italic.x = F) +
      ggtitle(gene)
    return(p)
  })
  if (length(pl) == 1) {
    p <- pl[[1]]
    ggsave(plot.out, p, device = cairo_pdf, height = sub.height, width = sub.width, dpi = 300)
  } else {
    p <- scFanc::wrap.plots.fanc(plot.list = pl, n.col = n.col,
                                 sub.width = sub.width, sub.height = sub.height,
                                 plot.out = plot.out)
  }
  invisible(p)
}

barkbase.rna <- function(file) {
  df <- read.table(file, header = T,
                   sep = "\t")
  df <- df %>% filter(GeneSymbol != "")
  df$GeneSymbol <- make.unique(df$GeneSymbol)
  rownames(df) <- df$GeneSymbol
  mat <- df[, -c(1:3)] %>% as.matrix()
  return(mat)
}

merge.bed.color <- function(beds, colors = NULL, out.rgbPeak, genome, ignore.strand = F,
                            bedToBigBed = "/opt/apps/kentUCSC/334/bedToBigBed") {
  # merge bed files into a rgbPeak formatted file.
  # the forth column and strand will be maintained.
  if (is.null(colors)) {
    colors <- RColorBrewer::colorRampPalette(c("red", "green"))(length(beds))
  }
  if (length(beds) != length(colors)) {
    stop("length(beds) != length(colors)")
  }
  colors <- sapply(colors, function(color) {
    if (!grepl(",", color)) {
      color <- col2rgb(color)[, 1] %>% paste0(collapse = ",")
    }
    return(color)
  })
  df.rgb <- lapply(seq_along(beds), function(i) {
    bed <- beds[i]
    if (length(readLines(bed)) == 0) {
      # return()
      bed <- data.frame(chr = "chr1", start = 1, end = 2)
    } else {
      bed <- utilsFanc::import.bed.fanc(bed)
    }
    bed <- bed[, 1:min(6, ncol(bed))]
    for (col in c("forth", "fifth")) {
      if (!col %in% colnames(bed)) {
        # bed[, col] <- paste0(rep("A", 10000), collapse = "")
        bed[, col] <- " " 
        # bed[, col] <- "A"
        # this is an invisible character copied from https://www.editpad.org/tool/invisible-character
      }
    }
    
    if (is.null(bed$strand) || ignore.strand) {
      bed$strand <- "."
    }
    
    bed$fifth <- 1
    bed$thickStart <- bed$start + 1
    bed$thickEnd <- bed$end
    bed$color <- colors[i]
    if (ncol(bed) != 9) {
      stop("ncol(bed) != 9")
    }
    return(bed)
  }) %>% do.call(rbind, .)
  dir.create(dirname(out.rgbPeak), showWarnings = F, recursive = T )
  df.rgb <- df.rgb %>% arrange(chr, start)
  temp.bed <- paste0(out.rgbPeak, ".bed")
  write.table(df.rgb, temp.bed, sep = "\t", quote = F, col.names = F, row.names = F)
  chrom.sizes <- paste0("~/genomes/", genome, "/", genome, ".chrom.sizes")
  cmd <- paste0(bedToBigBed, " ", temp.bed, " ", chrom.sizes, " ", out.rgbPeak)
  system(cmd)
  cat(utilsFanc::bash2ftp(filename = out.rgbPeak))
  cat("\n")
  return(out.rgbPeak)
}
