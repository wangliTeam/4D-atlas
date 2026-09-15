# ============================================================================
#TOME analysis
# Ref
#https://github.com/ChengxiangQiu/tome_code/blob/main/Section6_keyTF_Step1_increase_decrease.R
#https://github.com/ChengxiangQiu/tome_code/blob/main/Section6_keyTF_Step2_summarize_results.R


suppressPackageStartupMessages({
  library(Seurat)
  library(schard)   
  library(dplyr)
  library(harmony) 
  library(FNN)      
  library(networkD3) 
  library(Matrix)
  library(gtools)   
})

# ==============================================================================
# 1. Global configuration
# ==============================================================================

INPUT_ROOT  <- "../3DModel_Data"
OUTPUT_ROOT <- "../TOME"
OVERWRITE   <- TRUE  
#TF link
TF_FILES <- list(
  drosophila = "../Drosophila_melanogaster_TF.txt",
  mouse      = "../Mus_musculus_TF.txt",
  human      = "../Homo_sapiens_TF.txt"
)

aggregate_single_annotation <- function(file_list, target_col) {
  meta_states <- list()
  
  for (f in file_list) {
    tp <- sub("_alignment.h5ad$", "", basename(f))
    message(paste0("     -> Aggregating: ", tp, " ..."))
    
    # 1. Read
    obj <- tryCatch(schard::h5ad2seurat(f, use.raw = FALSE), error=function(e) NULL)
    if(is.null(obj)) { warning(paste("Read error:", f)); next }
    
    # 2. Check the columns
    if (!target_col %in% colnames(obj@meta.data)) {
      warning(paste("        [Skip] Column", target_col, "not found in", tp))
      next
    }
    
    # 3. Extract the matrix
    mat <- tryCatch(GetAssayData(obj, assay="RNA", layer="data"), error=function(e) GetAssayData(obj, assay="RNA", slot="data"))
    if (max(mat) < 50) mat <- expm1(mat)
    
    groups <- obj@meta.data[[target_col]]
    unique_groups <- unique(as.character(groups))
    unique_groups <- unique_groups[!is.na(unique_groups)]
    if (length(unique_groups) < 2) {next
    }
    avg_mat <- matrix(0, nrow = nrow(mat), ncol = length(unique_groups))
    colnames(avg_mat) <- unique_groups; rownames(avg_mat) <- rownames(mat)
    
    for (g in unique_groups) {
      cells <- which(groups == g)
      if (length(cells) > 1) avg_mat[, g] <- Matrix::rowMeans(mat[, cells, drop=FALSE])
      else avg_mat[, g] <- mat[, cells]
    }
    avg_mat <- log1p(avg_mat)
    i_j_x <- which(avg_mat != 0, arr.ind = TRUE)
    if(length(i_j_x) == 0) {
      avg_mat_sparse <- Matrix::Matrix(0, nrow=nrow(avg_mat), ncol=ncol(avg_mat), sparse=TRUE)
    } else {
      avg_mat_sparse <- Matrix::sparseMatrix(
        i = i_j_x[, 1],
        j = i_j_x[, 2],
        x = avg_mat[i_j_x],
        dims = c(nrow(avg_mat), ncol(avg_mat))
      )
    }
    dimnames(avg_mat_sparse) <- list(rownames(avg_mat), colnames(avg_mat))
    meta_df <- data.frame(State = colnames(avg_mat_sparse), Timepoint = tp, row.names = colnames(avg_mat_sparse))
    
    state_obj <- tryCatch({
      CreateSeuratObject(counts = avg_mat_sparse, meta.data = meta_df)
    }, error = function(e) {
      assay_data <- CreateAssayObject(counts = avg_mat_sparse)
      SeuratObject(assay = assay_data, meta.data = meta_df)
    })
    if (inherits(state_obj[["RNA"]], "Assay5")) {
      state_obj[["RNA"]]$data <- avg_mat_sparse
    } else {
      state_obj <- SetAssayData(state_obj, slot = "data", new.data = avg_mat_sparse)
    }
    
    state_obj <- FindVariableFeatures(state_obj, nfeatures = 2000, verbose = FALSE)
    state_obj <- ScaleData(state_obj, verbose = FALSE)
    meta_states[[tp]] <- state_obj
  }
  return(meta_states)
}

# --- A2. Trajectory calculation ---
calculate_transition <- function(obj_pre, obj_nex, k = 5) {
  merged <- merge(obj_pre, y = obj_nex, add.cell.ids = c("Pre", "Nex"))
  VariableFeatures(merged) <- intersect(VariableFeatures(obj_pre), VariableFeatures(obj_nex))
  merged <- ScaleData(merged, verbose = FALSE)
  
  n_pcs <- min(30, ncol(merged) - 1)
  if (n_pcs < 2) {
    emb <- t(GetAssayData(merged, slot = "scale.data"))
  } else {
    merged <- RunPCA(merged, npcs = n_pcs, verbose = FALSE)
    emb <- tryCatch({
      merged <- RunHarmony(merged, group.by.vars = "Timepoint", dims.use = 1:n_pcs, verbose = FALSE)
      Embeddings(merged, "harmony")
    }, error = function(e) Embeddings(merged, "pca"))
  }
  
  cells_pre <- colnames(merged)[merged$Timepoint == unique(obj_pre$Timepoint)]
  cells_nex <- colnames(merged)[merged$Timepoint == unique(obj_nex$Timepoint)]
  
  k_use <- min(k, length(cells_pre))
  if(k_use < 1) k_use <- 1
  
  knn_res <- get.knnx(data = emb[cells_pre, , drop=FALSE], query = emb[cells_nex, , drop=FALSE], k = k_use)
  
  trans_df <- data.frame()
  for (i in 1:nrow(knn_res$nn.index)) {
    w <- 1 / (knn_res$nn.dist[i, ] + 1e-6); w <- w / sum(w)
    tmp <- data.frame(
      Pre_State = merged@meta.data[cells_pre[knn_res$nn.index[i, ]], "State"],
      Next_State = merged@meta.data[cells_nex[i], "State"],
      Probability = w,
      Time_Start = unique(obj_pre$Timepoint),
      Time_End = unique(obj_nex$Timepoint)
    )
    trans_df <- rbind(trans_df, tmp)
  }
  
  if (nrow(trans_df) > 0) {
    final_df <- trans_df %>%
      group_by(Pre_State, Next_State, Time_Start, Time_End) %>%
      summarise(Sum_Prob = sum(Probability), .groups = 'drop') %>%
      group_by(Next_State) %>%
      slice_max(order_by = Sum_Prob, n = 1, with_ties = FALSE) %>%
      mutate(Probability = 1) %>% ungroup() %>% select(-Sum_Prob)
    return(final_df)
  } else { return(data.frame()) }
}

run_tome_tf_analysis <- function(transitions, obj_list, tf_list, direction = "increase") {
  results_list <- list()
  valid_tps <- names(obj_list)
  THRESH <- 0.25
  
  for (t_idx in 1:(length(valid_tps)-1)) {
    tp_curr <- valid_tps[t_idx]; tp_next <- valid_tps[t_idx+1]
    
    sub_edges <- transitions %>% filter(Time_Start == tp_curr, Time_End == tp_next, Probability > 0.5)
    if(nrow(sub_edges) == 0) next
    mat_pre <- tryCatch(GetAssayData(obj_list[[tp_curr]], layer="data"), error=function(e) GetAssayData(obj_list[[tp_curr]], slot="data"))
    mat_nex <- tryCatch(GetAssayData(obj_list[[tp_next]], layer="data"), error=function(e) GetAssayData(obj_list[[tp_next]], slot="data"))
    
    common_genes <- intersect(rownames(mat_pre), rownames(mat_nex))
    valid_tfs_here <- intersect(tf_list, common_genes)
    if(length(valid_tfs_here) < 2) next
    
    pre_states <- unique(sub_edges$Pre_State)
    for (parent in pre_states) {
      if (!parent %in% colnames(mat_pre)) next
      edge_sub <- sub_edges[sub_edges$Pre_State == parent, ]
      
      expr_parent <- mat_pre[valid_tfs_here, parent]
      for (j in 1:nrow(edge_sub)) {
        child <- edge_sub$Next_State[j]
        if (!child %in% colnames(mat_nex)) next
        
        expr_child <- mat_nex[valid_tfs_here, child]
        logfc_1 <- expr_child - expr_parent
        pass_1 <- if(direction == "increase") logfc_1 > THRESH else logfc_1 < -THRESH
        candidates <- names(logfc_1)[pass_1]
        if (length(candidates) == 0) next
        
        siblings <- intersect(edge_sub$Next_State[edge_sub$Next_State != child], colnames(mat_nex))
        logfc_2_vals <- rep(NA, length(candidates)); names(logfc_2_vals) <- candidates
        final_genes <- candidates
        
        if (length(siblings) > 0) {
          if (length(siblings) > 1) expr_sibs <- rowMeans(mat_nex[candidates, siblings, drop=FALSE])
          else expr_sibs <- mat_nex[candidates, siblings]
          logfc_2 <- expr_child[candidates] - expr_sibs
          keep <- if(direction == "increase") logfc_2 > THRESH else logfc_2 < -THRESH
          final_genes <- candidates[keep]
          logfc_2_vals <- logfc_2[keep]
        }
        if (length(final_genes) == 0) next
        
        results_list[[length(results_list) + 1]] <- data.frame(
          celltype = child, emergence_time = tp_next, key_TF = final_genes,
          edge = paste0(tp_curr, ":", parent, " -> ", tp_next, ":", child),
          Pre_State = parent, Next_State = child, Time_Start = tp_curr, Time_End = tp_next,
          avg_logFC_1 = logfc_1[final_genes], avg_logFC_2 = logfc_2_vals[final_genes],
          direction = direction
        )
      }
    }
  }
  res_df <- do.call(rbind, results_list)
  
  if (!is.null(res_df)) {
    res_df <- res_df %>% group_by(edge) %>%
      mutate(
        score_raw_1 = scale(abs(avg_logFC_1))[,1],
        score_raw_2 = ifelse(is.na(avg_logFC_2), NA, scale(abs(avg_logFC_2))[,1]),
        determined_score = ifelse(is.na(score_raw_2), score_raw_1, (score_raw_1 + score_raw_2)/2)
      ) %>%
      mutate(determined_score = ifelse(is.na(determined_score), abs(avg_logFC_1), determined_score)) %>%
      ungroup() %>% arrange(desc(determined_score))
  }
  return(res_df)
}
subfolders <- list.dirs(INPUT_ROOT, full.names = TRUE, recursive = FALSE)
if (length(subfolders) == 0) stop("No data folders found!")
subfolders <- subfolders[basename(subfolders) %in% target_datasets]

for (folder_path in subfolders) {
  dataset_name <- basename(folder_path)
  
  h5ad_files <- list.files(folder_path, pattern = "_alignment.h5ad$", full.names = TRUE)
  h5ad_files <- gtools::mixedsort(h5ad_files)
  names(h5ad_files) <- gsub("_alignment.h5ad$", "", basename(h5ad_files))
  first_obj <- tryCatch(schard::h5ad2seurat(h5ad_files[1], use.raw = FALSE), error=function(e) NULL)
  available_cols <- intersect(TARGET_ANNOTATIONS, colnames(first_obj@meta.data))

  for (current_col in available_cols) {   
    save_dir <- file.path(OUTPUT_ROOT, dataset_name)
    if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
    setwd(save_dir)
    file_suffix <- paste0("_", current_col)
    rds_name <- paste0("Aggregated_List", file_suffix, ".rds")
    if (file.exists(rds_name) && OVERWRITE == FALSE) {
      aggregated_list <- readRDS(rds_name)
    } else {
      aggregated_list <- aggregate_single_annotation(h5ad_files, current_col)
      saveRDS(aggregated_list, rds_name)
    }
    valid_tps <- names(aggregated_list)
    all_transitions <- data.frame()
    if (length(valid_tps) > 1) {
      for (i in 1:(length(valid_tps)-1)) {
        edges <- calculate_transition(aggregated_list[[valid_tps[i]]], aggregated_list[[valid_tps[i+1]]])
        all_transitions <- rbind(all_transitions, edges)
      }
    }
    write.csv(all_transitions, paste0("TOME_Trajectory_Results", file_suffix, ".csv"), row.names = FALSE)
    if(nrow(all_transitions) > 0) {
      nodes <- data.frame(name = unique(c(paste(all_transitions$Time_Start, all_transitions$Pre_State, sep=":"), 
                                          paste(all_transitions$Time_End, all_transitions$Next_State, sep=":"))))
      links <- all_transitions
      links$IDsource <- match(paste(links$Time_Start, links$Pre_State, sep=":"), nodes$name) - 1
      links$IDtarget <- match(paste(links$Time_End, links$Next_State, sep=":"), nodes$name) - 1
      sankey <- sankeyNetwork(Links = links, Nodes = nodes, Source = "IDsource", Target = "IDtarget", 
                              Value = "Probability", NodeID = "name", sinksRight = FALSE, fontSize = 12, nodeWidth = 30)
      saveNetwork(sankey, paste0("TOME_Tree_Sankey", file_suffix, ".html"), selfcontained = FALSE)
    }
    selected_tf_file <- NULL
    for (key in names(SPECIES_KEY_MAP)) {
      if (grepl(key, dataset_name, ignore.case = TRUE)) {
        selected_tf_file <- TF_FILES[[SPECIES_KEY_MAP[[key]]]]
        break
      }
    }
    
    if (!is.null(selected_tf_file) && file.exists(selected_tf_file)) {
      tryCatch({
        tf_data <- read.table(selected_tf_file, header=TRUE, as.is=TRUE, sep="\t")
        if("Symbol" %in% colnames(tf_data)) {
          species_tf <- as.vector(unique(tf_data$Symbol))
        } else {
          tf_data <- read.table(selected_tf_file, header=FALSE, as.is=TRUE, sep="\t")
          species_tf <- as.vector(unique(tf_data$V1))
        }
      }, error = function(e) {
        tf_data <- read.table(selected_tf_file, header=FALSE, as.is=TRUE, sep="\t")
        species_tf <- as.vector(unique(tf_data$V1))
      })
      
      if(length(species_tf) > 0) {
        df_inc <- run_tome_tf_analysis(all_transitions, aggregated_list, species_tf, "increase")
        df_dec <- run_tome_tf_analysis(all_transitions, aggregated_list, species_tf, "decrease")
        final_output <- rbind(df_inc, df_dec)
        if(!is.null(final_output)) {
          write.csv(final_output, paste0("TOME_KeyTF_Final", file_suffix, ".csv"), row.names = FALSE)
        }
      }
    }
  }
}