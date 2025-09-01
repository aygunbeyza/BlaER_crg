#!/usr/bin/env Rscript
# =========================
# Seurat + DoubletFinder batch pipeline
# =========================

suppressPackageStartupMessages({
  library(Seurat)
  library(DoubletFinder)
  library(ggplot2)
  library(dplyr)
  library(readr)
})

# ---------- User configuration (EDIT HERE) ----------
datasets <- list(
  list(name = "t0",
       in_dir  = "/t0_with_genefull_Solo.out/Gene/filtered",
       out_dir = "/t0_seurat_results"),
  list(name = "t120",
       in_dir  = "/t120_with_genefull_Solo.out/Gene/filtered",
       out_dir = "/t120_seurat_results")
  list(name = "t0_genefull",
       in_dir  = "/t0_with_genefull_Solo.out/GeneFull/filtered",
       out_dir = "/t0_genefull_seurat_results"),
  list(name = "t120_genefull",
       in_dir  = "/t120_with_genefull_Solo.out/GeneFull/filtered",
       out_dir = "/result_seurat/t120_genefull_seurat_results")
)

# Global analysis params
min_features   <- 200
max_features   <- 5000
max_percent_mt <- 20
n_pcs          <- 15
resolution     <- 0.3
doublet_rate   <- 0.085
set.seed(12345)

# ---------------------------------------------------

safe_percent_mt <- function(obj) {
  feats <- rownames(obj)
  pat <- if(any(grepl("^MT-", feats))) "^MT-" else "^mt-"
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = pat)
  obj
}

save_plot <- function(plot, file, w=7, h=6, dpi=300) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  ggsave(file, plot=plot, width=w, height=h, dpi=dpi)
}

process_one <- function(name, in_dir, out_dir) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  # per-dataset log
  log_file <- file.path(out_dir, paste0(name, "_pipeline_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
  zz <- file(log_file, open="wt")
  sink(zz, split=TRUE); sink(zz, type="message")
  on.exit({
    cat("\n--- SESSION INFO (END) ---\n"); print(utils::sessionInfo())
    sink(type="message"); sink(); close(zz)
    cat(sprintf("\nLog saved: %s\n", log_file))
  }, add=TRUE)

  cat(sprintf("=== [%s] START ===\nInput: %s\nOutput: %s\n", name, in_dir, out_dir))

  # Load 10x + QC
  counts <- Read10X(data.dir = in_dir)
  obj <- CreateSeuratObject(counts = counts, project = name, min.cells = 3, min.features = min_features)
  obj <- safe_percent_mt(obj)
  obj_filt <- subset(obj, subset = nFeature_RNA > min_features &
                               nFeature_RNA < max_features &
                               percent.mt < max_percent_mt)

  # SCT + PCA/UMAP/Neighbors/Clusters
  obj_filt <- SCTransform(obj_filt, verbose = FALSE)
  obj_filt <- RunPCA(obj_filt, assay = "SCT", npcs = max(30, n_pcs))
  obj_filt <- RunUMAP(obj_filt, dims = 1:n_pcs)
  obj_filt <- FindNeighbors(obj_filt, dims = 1:n_pcs)
  obj_filt <- FindClusters(obj_filt, resolution = resolution)
  DefaultAssay(obj_filt) <- "SCT"

  top10 <- head(VariableFeatures(obj_filt), 10)
  p_varL <- LabelPoints(plot = VariableFeaturePlot(obj_filt), points = top10, repel = TRUE) +
            ggtitle(sprintf("Variable features (Top 10) - %s", name))
  p_umap0 <- DimPlot(obj_filt, reduction="umap", label=TRUE) +
             ggtitle(sprintf("UMAP - filtering - %s", name))
  save_plot(p_varL,  file.path(out_dir, paste0(name, "_variable_features.png")), w=8, h=6)
  save_plot(p_umap0, file.path(out_dir, paste0(name, "_umap_filtering.png")))

  # DoubletFinder sweep
  obj_filt <- ScaleData(obj_filt, assay="SCT")
  sweep.res   <- paramSweep(obj_filt, PCs = 1:n_pcs, sct = TRUE)
  sweep.stats <- summarizeSweep(sweep.res, GT=FALSE)
  bcmvn       <- find.pK(sweep.stats)

  write_csv(bcmvn, file.path(out_dir, paste0(name, "_bcmvn.csv")))
  p_bc <- ggplot(bcmvn, aes(pK, BCmetric, group=1)) + geom_point() + geom_line() +
          ggtitle(sprintf("DoubletFinder - pK vs BCmetric (%s)", name))
  save_plot(p_bc, file.path(out_dir, paste0(name, "_doubletfinder_pK_BCmetric.png")), w=7, h=5)

  bcmvn$pK <- as.numeric(as.character(bcmvn$pK))
  pK_value <- bcmvn$pK[ which.max(bcmvn$BCmetric) ]
  cat(sprintf("Chosen pK_value = %.4f\n", pK_value))

  annotations    <- obj_filt$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)
  nExp_poi       <- round(doublet_rate * nrow(obj_filt@meta.data))
  nExp_poi.adj   <- round(nExp_poi * (1 - homotypic.prop))
  cat(sprintf("nExp_poi = %d; nExp_poi.adj = %d; homotypic.prop = %.3f\n",
              nExp_poi, nExp_poi.adj, homotypic.prop))

  # DoubletFinder classify
  obj_df <- doubletFinder(
    obj_filt, PCs=1:n_pcs, pN=0.25, pK=pK_value, nExp=nExp_poi.adj, reuse.pANN=NULL, sct=TRUE
  )
  df_col <- grep("^DF\\.classifications", colnames(obj_df@meta.data), value=TRUE)[1]
  if (is.na(df_col)) stop("No DF.classifications column found.")

  p_umap_df <- DimPlot(obj_df, reduction="umap", group.by=df_col) +
               ggtitle(sprintf("DoubletFinder - UMAP classification (%s; pK=%.3f)", name, pK_value))
  save_plot(p_umap_df, file.path(out_dir, paste0(name, "_doubletfinder_umap_classification.png")))

  df_counts <- as.data.frame(table(obj_df@meta.data[[df_col]]))
  write_csv(df_counts, file.path(out_dir, paste0(name, "_doubletfinder_class_counts.csv")))
  p_bar <- ggplot(df_counts, aes(x=Var1, y=Freq, fill=Var1)) +
           geom_bar(stat="identity") + theme_minimal() +
           xlab("Classification") + ylab("Number of cells") +
           ggtitle(sprintf("Doublet vs Singlet counts - %s", name))
  save_plot(p_bar, file.path(out_dir, paste0(name, "_doubletfinder_counts_barplot.png")), w=6, h=5)

  # Subset Singlets & re-run
  obj_sing <- subset(obj_df, subset = !!as.name(df_col) == "Singlet")
  DefaultAssay(obj_sing) <- "SCT"
  obj_sing <- RunPCA(obj_sing, assay="SCT", npcs = max(30, n_pcs))
  obj_sing <- RunUMAP(obj_sing, dims = 1:n_pcs)
  obj_sing <- FindNeighbors(obj_sing, dims = 1:n_pcs)
  obj_sing <- FindClusters(obj_sing, resolution = resolution)

  top10_after   <- head(VariableFeatures(obj_sing), 10)
  p_var_afterL  <- LabelPoints(plot=VariableFeaturePlot(obj_sing), points=top10_after, repel=TRUE) +
                   ggtitle(sprintf("Variable features (Top 10) - %s (after DoubletFinder)", name))
  p_umap_after  <- DimPlot(obj_sing, reduction="umap", label=TRUE) +
                   ggtitle(sprintf("UMAP - %s (after DoubletFinder)", name))
  save_plot(p_var_afterL, file.path(out_dir, paste0(name, "_doubletfinder_after_variable_features.png")), w=8, h=6)
  save_plot(p_umap_after, file.path(out_dir, paste0(name, "_doubletfinder_after_umap.png")))

# ---- Top10 marker heatmap + gene list saver ----
make_cluster_heatmap <- function(obj, out_dir, n_top = 10, group_by = "seurat_clusters",
                                 logfc_cut = 1) {
  DefaultAssay(obj) <- "SCT"

  # Küme kimlikleri
  if (!group_by %in% colnames(obj@meta.data)) {
    stop(sprintf("Column '%s' not found in meta.data", group_by))
  }
  Idents(obj) <- obj@meta.data[[group_by]]

  # Markerları bul (pozitif markerlar)
  markers <- FindAllMarkers(
    obj, only.pos = TRUE, logfc.threshold = 0, min.pct = 0.1, verbose = FALSE
  )

  # Bazı Seurat sürümlerinde sütun adı farklı olabilir
  if (!"avg_log2FC" %in% colnames(markers) && "avg_logFC" %in% colnames(markers)) {
    markers$avg_log2FC <- markers$avg_logFC
  }

  # Filtrele ve cluster bazında top-N seç
  top_markers <- markers |>
    dplyr::filter(avg_log2FC > logfc_cut) |>
    dplyr::group_by(cluster) |>
    dplyr::slice_head(n = n_top) |>
    dplyr::ungroup()

  # Isı haritası (gen isimleri unique olmalı)
  feats <- unique(top_markers$gene)
  p_heat <- DoHeatmap(obj, features = feats, group.by = group_by) + NoLegend()

  # Kayıtlar
  heat_f   <- file.path(out_dir, sprintf("%s_top%d_markers_heatmap.png", name, n_top))
  table_f  <- file.path(out_dir, sprintf("%s_top%d_markers_per_cluster.csv", name, n_top))
  save_plot(p_heat, heat_f, w = 10, h = 8)

  # Top marker tabloyu kaydet
  readr::write_csv(top_markers, table_f)

  # Cluster başına .txt (genler alt alta)
  clus_dir <- file.path(out_dir, sprintf("%s_top%d_marker_lists", name, n_top))
  dir.create(clus_dir, recursive = TRUE, showWarnings = FALSE)
  splitted <- split(top_markers$gene, top_markers$cluster)
  for (cl in names(splitted)) {
    fn <- file.path(clus_dir, sprintf("cluster_%s_top%d_genes.txt", cl, n_top))
    writeLines(splitted[[cl]], fn)
  }

  invisible(list(heatmap_file = heat_f, table_file = table_f, dir = clus_dir))
}

# ---- Top10 marker heatmap + listeler (DF sonrası singlet obje) ----
make_cluster_heatmap(obj_sing, out_dir, n_top = 10, group_by = "seurat_clusters", logfc_cut = 1)


  # Save RDS
  saveRDS(obj_filt, file.path(out_dir, paste0(name, "_filtered_sct.rds")))
  saveRDS(obj_df,   file.path(out_dir, paste0(name, "_df_annotated.rds")))
  saveRDS(obj_sing, file.path(out_dir, paste0(name, "_singlet_after_doubletfinder.rds")))

  cat(sprintf("=== [%s] DONE ===\n", name))
}

# Run all
for (ds in datasets) {
  tryCatch({
    process_one(name=ds$name, in_dir=ds$in_dir, out_dir=ds$out_dir)
  }, error=function(e){
    message(sprintf("\n[ERROR] Dataset '%s' failed: %s", ds$name, conditionMessage(e)))
  })
}
