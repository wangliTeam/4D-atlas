# ============================================================================
#Temporal expression clustering (Mfuzz)
#ref
#https://github.com/BGI-DEV-REG/ARTISTA/blob/main/code/Figure5I_J_Mfuzz.R
suppressPackageStartupMessages({
  library(Seurat)
  library(schard)   
  library(dplyr)
  library(Matrix)
  library(gtools)   
  library(ggplot2)
  library(reshape2)
  library(Biobase)
  library(Mfuzz)
})

CLUSTER_NUM     <- 12  
MIN_TIMEPOINTS  <- 3  
MIN_EXPR_THRESH <- 0.1 
TOP_VAR_GENES   <- 3000 

COLOR_BG_LINE   <- "#7fbf7b" 
COLOR_CENTER    <- "#F7941D" 

run_mfuzz_analysis <- function(expr_mat, out_dir, prefix) {

  expr_mat <- expr_mat[rowMeans(expr_mat) > MIN_EXPR_THRESH, , drop=FALSE]
  if (nrow(expr_mat) < 50) return(NULL) 

  if (!is.null(TOP_VAR_GENES) && nrow(expr_mat) > TOP_VAR_GENES) {
    gene_sd <- apply(expr_mat, 1, sd)
    expr_mat <- expr_mat[order(gene_sd, decreasing = TRUE)[1:TOP_VAR_GENES], ]
  }
  eset <- ExpressionSet(assayData = as.matrix(expr_mat))
  eset_s <- standardise(eset)
  m_val <- tryCatch(mestimate(eset_s), error=function(e) 1.5)
  set.seed(123) 
  cl <- mfuzz(eset_s, c = CLUSTER_NUM, m = m_val)
  centers <- as.data.frame(cl$centers)
  centers$Cluster <- paste0("Cluster ", 1:nrow(centers))
  centers$Cluster <- factor(centers$Cluster, levels = paste0("Cluster ", 1:nrow(centers)))
  centers_long <- melt(centers, id.vars = "Cluster", variable.name = "Timepoint", value.name = "Expression")
  mem_df <- as.data.frame(cl$membership)
  mem_df$Gene <- rownames(mem_df)
  mem_df$Max_Mem <- apply(mem_df[, 1:CLUSTER_NUM], 1, max)
  mem_df$Best_Cluster <- apply(mem_df[, 1:CLUSTER_NUM], 1, which.max)
  core_genes <- mem_df$Gene[mem_df$Max_Mem > 0.3]
  expr_std <- as.data.frame(exprs(eset_s))
  expr_std$Gene <- rownames(expr_std)
  expr_long <- melt(expr_std[core_genes, ], id.vars = "Gene", variable.name = "Timepoint", value.name = "Expression")
  cluster_info <- data.frame(Gene = mem_df$Gene, Cluster = paste0("Cluster ", mem_df$Best_Cluster))
  expr_long <- merge(expr_long, cluster_info, by = "Gene")
  expr_long$Cluster <- factor(expr_long$Cluster, levels = levels(centers$Cluster))
  
  p <- ggplot() +
    geom_line(data = plot_bg, aes(x = Timepoint, y = Expression, group = Gene), 
              color = COLOR_BG_LINE, size = 0.1, alpha = 0.4) +
    geom_line(data = centers_long, aes(x = Timepoint, y = Expression, group = Cluster), 
              color = COLOR_CENTER, size = 1.2) +
    geom_point(data = centers_long, aes(x = Timepoint, y = Expression, group = Cluster), 
               color = COLOR_CENTER, size = 1.5) +
    
    facet_wrap(~Cluster, ncol = 3, scales = "free_y") + 
    labs(title = prefix, x = "", y = "Standardized Expression") +
    my_theme + axis_opts
  h_val <- if(CLUSTER_NUM > 8) 12 else 10
  ggsave(file.path(out_dir, paste0(prefix, "_Trend_Plot.pdf")), p, width = 10, height = h_val)

  write.csv(mem_df, file.path(out_dir, paste0(prefix, "_Gene_Clusters.csv")), row.names = FALSE)

}
subfolders <- list.dirs(INPUT_ROOT, full.names = TRUE, recursive = FALSE)
subfolders <- subfolders[basename(subfolders) %in% target_datasets]

for (folder_path in subfolders) {
  dataset_name <- basename(folder_path)
  h5ad_files <- list.files(folder_path, pattern = "_alignment.h5ad$", full.names = TRUE)
  h5ad_files <- gtools::mixedsort(h5ad_files) 
  timepoints <- gsub("_alignment.h5ad$", "", basename(h5ad_files))
  names(h5ad_files) <- timepoints
  global_means_list <- list()
  anno_means_list <- list(); for(c in TARGET_COLS) anno_means_list[[c]] <- list()
  all_genes <- NULL
  for (tp in timepoints) {
    obj <- tryCatch(schard::h5ad2seurat(h5ad_files[[tp]], use.raw = FALSE), error=function(e) NULL)
    if(is.null(obj)) next
    mat <- tryCatch(GetAssayData(obj, assay="RNA", layer="data"), error=function(e) GetAssayData(obj, assay="RNA", slot="data"))
    if (max(mat) < 50) mat <- expm1(mat)
    if (is.null(all_genes)) all_genes <- rownames(mat) else all_genes <- intersect(all_genes, rownames(mat))
    global_means_list[[tp]] <- Matrix::rowMeans(mat)
    meta <- obj@meta.data
    for (col in TARGET_COLS) {
      if (col %in% colnames(meta)) {
        groups <- as.character(meta[[col]])
        unique_gps <- unique(groups[!is.na(groups)])
        tmp_mat <- matrix(NA, nrow = nrow(mat), ncol = length(unique_gps))
        colnames(tmp_mat) <- unique_gps; rownames(tmp_mat) <- rownames(mat)
        for (g in unique_gps) {
          cells <- which(groups == g)
          if (length(cells) > 1) {
            tmp_mat[, g] <- Matrix::rowMeans(mat[, cells, drop=FALSE])
          } else {
            tmp_mat[, g] <- mat[, cells]
          }
        }
        anno_means_list[[col]][[tp]] <- tmp_mat
      }
    }
    rm(obj, mat); gc()
  }
  save_root <- file.path(OUTPUT_ROOT, dataset_name)
  if (!dir.exists(save_root)) dir.create(save_root, recursive = TRUE)
  global_mat <- do.call(cbind, lapply(global_means_list, function(x) x[all_genes]))
  global_mat <- log1p(global_mat) 
  run_mfuzz_analysis(global_mat, save_root, "Global_AllCells")
  for (col in TARGET_COLS) {
    tps_data <- anno_means_list[[col]]
    if (length(tps_data) < MIN_TIMEPOINTS) next
    col_save_dir <- file.path(save_root, col)
    all_groups <- unique(unlist(lapply(tps_data, colnames)))
    for (grp in all_groups) {
      valid_tps <- names(tps_data)[sapply(tps_data, function(x) grp %in% colnames(x))]
      
      if (length(valid_tps) >= MIN_TIMEPOINTS) {
        expr_list <- list()
        for (tp in valid_tps) expr_list[[tp]] <- tps_data[[tp]][all_genes, grp]
        
        grp_mat <- do.call(cbind, expr_list)
        grp_mat <- log1p(grp_mat) 
        
        safe_name <- gsub("[^a-zA-Z0-9_]", "_", grp)
        run_mfuzz_analysis(grp_mat, col_save_dir, safe_name)
      }
    }
  }
}


