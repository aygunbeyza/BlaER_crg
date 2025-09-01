suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
  library(ggplot2)
  library(data.table)
})

lr_t0   <- read.table("t0_gene_exp_matrix.tsv.gz", 
                      header = TRUE, sep = "\t")
lr_t120 <- read.table("/t120_gene_exp_matrix.tsv.gz", 
                      header = TRUE, sep = "\t")

to_pseudobulk_by_ensembl <- function(df) {
  meta_cols <- c("Geneid", "gene_id", "gene_name", "Chr", "Start", "End", "Strand", "Length")
  cnt_df <- df[, !(names(df) %in% meta_cols), drop = FALSE]

  # clean NA
  mat <- as.matrix(sapply(cnt_df, function(v) as.numeric(as.character(v))))
  mat[is.na(mat)] <- 0

  # delete version
  ens_ids <- gsub("\\..*", "", df$Geneid)

  # sun the same ensenbl
  agg <- rowsum(mat, group = ens_ids, reorder = FALSE)

  v <- rowSums(agg)
  names(v) <- rownames(agg)
  v
}

lr_t0_pb   <- to_pseudobulk_by_ensembl(lr_t0)
lr_t120_pb <- to_pseudobulk_by_ensembl(lr_t120)

apply_ensembl_rownames <- function(count_matrix, features_path) {
  # read features.tsv.gz 
  features <- read.delim(features_path, header = FALSE)
  colnames(features) <- c("gene_id", "gene_symbol", "feature_type")
  features$gene_id <- gsub("\\..*", "", features$gene_id)
  if (nrow(features) != nrow(count_matrix)) {
    warning("Satır sayısı uyuşmuyor: features.tsv.gz ile count matrisin farklı!")
  }
  # change the names
  old_rownames <- rownames(count_matrix)
  rownames(count_matrix) <- features$gene_id

  # find empty 
  na_rows <- is.na(rownames(count_matrix)) | rownames(count_matrix) == ""
  if (any(na_rows)) {
    cat("Uyarı: ", sum(na_rows), " satırın Ensembl ID'si boş veya NA idi, atılıyor.\n")
    count_matrix <- count_matrix[!na_rows, , drop = FALSE]
  }
  count_matrix <- rowsum(as.matrix(count_matrix), group = rownames(count_matrix), reorder = FALSE)
  return(count_matrix)
}


sr_t0_counts <- Read10X(data.dir = "/t0_with_genefull_Solo.out/Gene/filtered/")
sr_t0_counts <- apply_ensembl_rownames(sr_t0_counts, "/t0_with_genefull_Solo.out/Gene/filtered/features.tsv.gz")
sr_t0_genefull_counts <- Read10X("/t0_with_genefull_Solo.out/GeneFull/filtered/")
sr_t0_genefull_counts <- apply_ensembl_rownames(sr_t0_genefull_counts, "/t0_with_genefull_Solo.out/GeneFull/filtered/features.tsv.gz")
sr_t120_counts <- Read10X("/t120_with_genefull_Solo.out/Gene/filtered/")
sr_t120_counts <- apply_ensembl_rownames(sr_t120_counts, "/t120_with_genefull_Solo.out/Gene/filtered/features.tsv.gz")
sr_t120_genefull_counts <- Read10X("/t120_with_genefull_Solo.out/GeneFull/filtered/")
sr_t120_genefull_counts <- apply_ensembl_rownames(sr_t120_genefull_counts, "/t120_with_genefull_Solo.out/GeneFull/filtered/features.tsv.gz")

# pseudobulk vectors
sr_t0_pb <- Matrix::rowSums(sr_t0_counts)
sr_t0_genefull_pb <- Matrix::rowSums(sr_t0_genefull_counts)
sr_t120_pb <- Matrix::rowSums(sr_t120_counts)
sr_t120_genefull_pb <- Matrix::rowSums(sr_t120_genefull_counts)


isoquant_t120 <- read.csv("/result_isoquant_t120/OUT/OUT.gene_counts.tsv", sep = "\t")

isoquant_t0 <- read.csv("/result_isoquant_t0/OUT/OUT.gene_counts.tsv", sep = "\t")

star_t120_bulk <- data.table::fread("/gene_counts/t120_gene_counts.txt", sep = "\t")

star_t0_bulk <- data.table::fread("/gene_counts/t0_gene_counts.txt", sep = "\t")

# 5. normalization
to_cpm <- function(counts) {
  1e6 * counts / sum(counts)
}

to_log_cpm <- function(counts) {
  log2(to_cpm(counts) + 1)
}

# for long read
lr_t0_logcpm   <- to_log_cpm(lr_t0_pb)
lr_t120_logcpm <- to_log_cpm(lr_t120_pb)

# for short read
sr_t0_logcpm            <- to_log_cpm(sr_t0_pb)
sr_t0_genefull_logcpm   <- to_log_cpm(sr_t0_genefull_pb)
sr_t120_logcpm          <- to_log_cpm(sr_t120_pb)
sr_t120_genefull_logcpm <- to_log_cpm(sr_t120_genefull_pb)

# isoquant bulk
isoquant_t0_vec <- isoquant_t0[[2]]  # count sütunu
names(isoquant_t0_vec) <- gsub("\\..*", "", isoquant_t0[[1]])  # #feature_id sütunu
isoquant_t0_logcpm <- to_log_cpm(isoquant_t0_vec)

isoquant_t120_vec <- isoquant_t120[[2]]
names(isoquant_t120_vec) <- gsub("\\..*", "", isoquant_t120[[1]])
isoquant_t120_logcpm <- to_log_cpm(isoquant_t120_vec)

# star bulk
star_t0_bulk_vec <- star_t0_bulk[[ncol(star_t0_bulk)]]
names(star_t0_bulk_vec) <- gsub("\\..*", "", star_t0_bulk$Geneid)
star_t0_bulk_logcpm <- to_log_cpm(star_t0_bulk_vec)

star_t120_bulk_vec <- star_t120_bulk[[ncol(star_t120_bulk)]]
names(star_t120_bulk_vec) <- gsub("\\..*", "", star_t120_bulk$Geneid)
star_t120_bulk_logcpm <- to_log_cpm(star_t120_bulk_vec)

# 6. corelation plot function
library(dplyr)
library(ggplot2)

# all normalize logCPM vectors
samples <- list(
  "lr_t0"            = lr_t0_logcpm,
  "lr_t120"          = lr_t120_logcpm,
  "sr_t0"            = sr_t0_logcpm,
  "sr_t0_full"       = sr_t0_genefull_logcpm,
  "sr_t120"          = sr_t120_logcpm,
  "sr_t120_full"     = sr_t120_genefull_logcpm,
  "isoquant_t0"      = isoquant_t0_logcpm,
  "isoquant_t120"    = isoquant_t120_logcpm,
  "star_bulk_t0"     = star_t0_bulk_logcpm,
  "star_bulk_t120"   = star_t120_bulk_logcpm
)

# table for results
correlation_summary <- data.frame()

# corelation function
analyze_correlation <- function(x_vec, y_vec, x_name, y_name) {
  common_genes <- intersect(names(x_vec), names(y_vec))
  expressed_genes <- common_genes[(x_vec[common_genes] > 0) | (y_vec[common_genes] > 0)]

  x_vals <- x_vec[expressed_genes]
  y_vals <- y_vec[expressed_genes]

  cor_val <- cor(x_vals, y_vals, method = "pearson")
  p <- ggplot(data.frame(x = x_vals, y = y_vals), aes(x = x, y = y)) +
    geom_point(alpha = 0.5) +
    labs(
      title = paste("Pearson Corr:", round(cor_val, 3)),
      subtitle = paste0(x_name, " vs ", y_name),
      x = paste0(x_name, " (logCPM)"),
      y = paste0(y_name, " (logCPM)")
    ) +
    coord_fixed() +
    xlim(0, max(c(x_vals, y_vals))) +
    ylim(0, max(c(x_vals, y_vals))) +
    theme_bw()

  print(p)

  correlation_summary <<- rbind(
    correlation_summary,
    data.frame(
      sample_x = x_name,
      sample_y = y_name,
      used_genes = length(expressed_genes),
      pearson_correlation = round(cor_val, 4)
    )
  )
}

# for all combinations
sample_names <- names(samples)
for (i in 1:(length(sample_names) - 1)) {
  for (j in (i + 1):length(sample_names)) {
    analyze_correlation(
      x_vec = samples[[i]],
      y_vec = samples[[j]],
      x_name = sample_names[i],
      y_name = sample_names[j]
    )
  }
}

correlation_summary
write.table(
  correlation_summary,
  file = "correlation_summary.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
