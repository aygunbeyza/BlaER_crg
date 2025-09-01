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
  list(name="t0",           out_dir="/t0_seurat_results"),
  list(name="t120",         out_dir="/t120_seurat_results"),
  list(name="t0_genefull",  out_dir="/t0_genefull_seurat_results"),
  list(name="t120_genefull",out_dir="/t120_genefull_seurat_results")
)

# ---- biomaRt ----
mart <- useEnsembl("genes", dataset="hsapiens_gene_ensembl")

is_ensembl <- function(x) grepl("^ENSG\\d+", x)

map_ensg_to_symbol <- function(ensg) {
  ensg <- unique(ensg)
  if (length(ensg) == 0) return(character(0))
  tb <- getBM(
    attributes = c("ensembl_gene_id","hgnc_symbol","external_gene_name"),
    filters    = "ensembl_gene_id",
    values     = ensg,
    mart       = mart
  )
  tb$symbol <- ifelse(nzchar(tb$hgnc_symbol), tb$hgnc_symbol, tb$external_gene_name)
  lut <- setNames(tb$symbol, tb$ensembl_gene_id)
  unname(lut[ensg])
}

# ---- filters (SYMBOL bazlı) ----
is_mt_symbol <- function(sym) grepl("^MT-", sym)
is_ribosomal_protein <- function(sym) grepl("^(RPS|RPL)[0-9A-Z]+$", sym)
is_rRNA_symbol <- function(sym) {
  # İnsan rRNA ailesi + mt-rRNA ve genel 'rRNA' yakalama
  grepl("^RNA(18|28|45)S$", sym) |
  grepl("^RNA5-?8S$", sym, ignore.case = TRUE) |
  grepl("^RNA5S$",  sym, ignore.case = TRUE) |
  grepl("^MT-RNR",  sym, ignore.case = TRUE) |
  grepl("rRNA",     sym, ignore.case = TRUE)
}

for (ds in datasets) {
  rds_path <- file.path(ds$out_dir, paste0(ds$name, "_singlet_after_doubletfinder.rds"))
  message(sprintf("[%s] loading: %s", ds$name, rds_path))
  obj <- readRDS(rds_path)

  counts <- GetAssayData(obj, assay="RNA", slot="counts")
  g_all  <- rownames(counts)[Matrix::rowSums(counts > 0) > 0]

    # ENSG -> SYMBOL map
  ensg <- g_all[is_ensembl(g_all)]
  symb <- g_all[!is_ensembl(g_all)]
  mapped <- map_ensg_to_symbol(ensg)

  report <- data.frame(
    original = c(ensg, symb),
    symbol   = c(mapped, symb),
    status   = c(ifelse(!is.na(mapped) & nzchar(mapped), "mapped_from_ENSG", "unmapped_ENSG"),
                 rep("already_symbol", length(symb))),
    stringsAsFactors = FALSE
  )

  # ---- temiz sembol listesi (NA/boş at) ----
  symbols_all <- sort(unique(report$symbol[!is.na(report$symbol) & nzchar(report$symbol)]))

  # ---- filtreler: MT / RPS-RPL / rRNA çıkar ----
  to_drop <- symbols_all[
    is_mt_symbol(symbols_all) | is_ribosomal_protein(symbols_all) | is_rRNA_symbol(symbols_all)
  ]
  symbols_filtered <- sort(setdiff(symbols_all, to_drop))

  # ---- dosya yolları ----
  out_syms_all  <- file.path(ds$out_dir, paste0(ds$name, "_expressed_genes.symbols.txt"))
  out_rep	<- file.path(ds$out_dir, paste0(ds$name, "_mapping_report.tsv"))
  out_unmapped  <- file.path(ds$out_dir, paste0(ds$name, "_unmapped_ENSG.txt"))

  out_syms_filt <- file.path(ds$out_dir, paste0(ds$name, "_expressed_genes.symbols.filtered.txt"))
  out_rep_filt  <- file.path(ds$out_dir, paste0(ds$name, "_mapping_report.filtered.tsv"))
  out_dropped   <- file.path(ds$out_dir, paste0(ds$name, "_dropped_genes.MT_RP_rRNA.txt"))

  # ---- yazımlar (ham) ----
  writeLines(symbols_all, out_syms_all)
  write_tsv(report, out_rep)
  writeLines(sort(unique(report$original[report$status=="unmapped_ENSG"])), out_unmapped)

  # ---- yazımlar (filtrelenmiş) ----
  # filtrelenmiş rapor: yalnızca kalan semboller
  report_filt <- report %>%
    filter(symbol %in% symbols_filtered & !(status == "unmapped_ENSG"))
  writeLines(symbols_filtered, out_syms_filt)
  write_tsv(report_filt, out_rep_filt)
  # düşenleri de kaydet (ne çıkarıldı görmek için)
  writeLines(sort(unique(to_drop)), out_dropped)

  message(sprintf("[%s] symbols: %d (all) -> %d (filtered).", ds$name, length(symbols_all), length(symbols_filtered)))
  message(sprintf("[%s] wrote:\n  - %s\n  - %s\n  - %s\n  - %s\n  - %s\n  - %s",
                  ds$name, out_syms_all, out_rep, out_unmapped, out_syms_filt, out_rep_filt, out_dropped))
}

