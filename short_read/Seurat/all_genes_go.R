#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(readr)
  library(biomaRt)
})


# ---- dataset list ----
datasets <- list(
  list(name="t0",   out_dir="/users/rg/baygun/BlaER_crg/short_read/Seurat/result_seurat/t0_seurat_results"),
  list(name="t120", out_dir="/users/rg/baygun/BlaER_crg/short_read/Seurat/result_seurat/t120_seurat_results")
)

# biomaRt bağlan
mart <- useEnsembl("genes", dataset="hsapiens_gene_ensembl")

is_ensembl <- function(x) grepl("^ENSG\\d+$", x)

map_ensg_to_symbol <- function(ensg) {
  if (length(ensg) == 0) return(character(0))
  tb <- getBM(attributes=c("ensembl_gene_id","hgnc_symbol","external_gene_name"),
              filters="ensembl_gene_id", values=unique(ensg), mart=mart)
  tb$symbol <- ifelse(nzchar(tb$hgnc_symbol), tb$hgnc_symbol, tb$external_gene_name)
  lut <- setNames(tb$symbol, tb$ensembl_gene_id)
  unname(lut[ensg])
}

for (ds in datasets) {
  rds_path <- file.path(ds$out_dir, paste0(ds$name, "_singlet_after_doubletfinder.rds"))
  message(sprintf("[%s] loading: %s", ds$name, rds_path))
  obj <- readRDS(rds_path)

  counts <- GetAssayData(obj, assay="RNA", slot="counts")
  g_all  <- rownames(counts)[rowSums(counts) > 0]

  # ENSG'leri SYMBOL'a çevir
  ensg <- g_all[is_ensembl(g_all)]
  symb <- g_all[!is_ensembl(g_all)]
  mapped <- map_ensg_to_symbol(ensg)

  report <- tibble(
    original = c(ensg, symb),
    symbol   = c(mapped, symb),
    status   = c(ifelse(!is.na(mapped), "mapped_from_ENSG","unmapped_ENSG"),
                 rep("already_symbol", length(symb)))
  )

  # nihai symbol listesi
  symbols <- sort(unique(report$symbol[!is.na(report$symbol) & report$symbol != ""]))

  # dosyalar
  out_syms <- file.path(ds$out_dir, paste0(ds$name, "_expressed_genes.symbols.txt"))
  out_rep  <- file.path(ds$out_dir, paste0(ds$name, "_mapping_report.tsv"))
  out_unm  <- file.path(ds$out_dir, paste0(ds$name, "_unmapped_ENSG.txt"))

  writeLines(symbols, out_syms)
  write_tsv(report, out_rep)
  writeLines(sort(unique(report$original[report$status=="unmapped_ENSG"])), out_unm)

  message(sprintf("[%s] wrote %d symbols -> %s", ds$name, length(symbols), out_syms))
}
