#' @export
process_grn_analysis <- function(obj, correlation_threshold = 0.4) {
  # Convert the RNA assay to the appropriate class
  obj[["RNA"]] <- as(obj[["RNA"]], "Assay")
  
  # Select TFs and generate heatmap
  res_tfs <- SelectTFs(object = obj, return.heatmap = TRUE)
  df.cor <- res_tfs$tfs
  ht <- res_tfs$heatmap
  
  # Select genes and generate heatmap
  res_genes <- SelectGenes(object = obj, labelTop1 = 0, labelTop2 = 0)
  df.p2g <- res_genes$p2g
  ht2 <- res_genes$heatmap
  
  # Compute TF-gene correlations
  tf.gene.cor <- GetTFGeneCorrelation(
    object = obj,
    tf.use = df.cor$tfs,
    gene.use = unique(df.p2g$gene),
    tf.assay = "chromvar",
    gene.assay = "RNA",
    trajectory.name = "Trajectory"
  )
  
  # Generate heatmaps for TF-gene correlations
  ht3 <- GRNHeatmap(tf.gene.cor)
  ht4 <- GRNHeatmap(tf.gene.cor, tf.timepoint = df.cor$time_point)
  
  # Prepare motif matching matrix
  motif.matching <- obj@assays$ATAC@motifs@data
  colnames(motif.matching) <- obj@assays$ATAC@motifs@motif.names
  motif.matching <- motif.matching[unique(df.p2g$peak), unique(tf.gene.cor$tf)]
  
  # Generate the Gene Regulatory Network (GRN)
  df.grn <- GetGRN(
    motif.matching = motif.matching,
    df.cor = tf.gene.cor,
    df.p2g = df.p2g
  )
  
  # Filter GRN based on correlation threshold and select key columns
  df.grn2 <- df.grn %>%
    subset(correlation > correlation_threshold) %>%
    dplyr::select(c(tf, gene, correlation))
  
  # Return results as a list
  return(list(
    heatmap_tf = ht,
    heatmap_genes = ht2,
    heatmap_tf_gene = ht3,
    heatmap_tf_gene_timepoint = ht4,
    grn = df.grn,
    grn_filtered = df.grn2
  ))
}


#' @export
coembed_data_function <- function(ATAC_objs, RNA_objs, gene.activities) {
    # Convert assays to "Assay5" format for each object
    for (i in 1:length(ATAC_objs)) {
        if ("peaks" %in% names(ATAC_objs[[i]]@assays)) {
            ATAC_objs[[i]][["peaks"]] <- as(ATAC_objs[[i]][["peaks"]], "Assay5")
        } else {
            message(paste0("No 'peaks' assay found in ATAC object ", i))
        }
        if ("RNA" %in% names(RNA_objs[[i]]@assays)) {
            RNA_objs[[i]][["RNA"]] <- as(RNA_objs[[i]][["RNA"]], "Assay5")
        } else {
            message(paste0("No 'RNA' assay found in RNA object ", i))
        }
    }
    
    # Initialize lists to store results
    co_embed_list <- list()
    pair_list <- list()
    final_pair_list <- list()
    
    # Co-embedding and harmonization loop
    for (i in 1:length(ATAC_objs)) {
        message(paste0("Co-embedding object ", i))
        
        # Filter gene activities based on shared column names with the ATAC object
        gene.activities_f <- gene.activities[, colnames(gene.activities) %in% colnames(ATAC_objs[[i]])]
        
        # Set Default Assay for both ATAC and RNA objects
        DefaultAssay(ATAC_objs[[i]]) <- 'RNA'
        DefaultAssay(RNA_objs[[i]]) <- 'RNA'
        
        # Perform co-embedding
        co_embed_list[[i]] <- CoembedData(
            RNA_objs[[i]],
            ATAC_objs[[i]],
            gene.activities_f,
            reference.assay = "RNA",
            weight.reduction = "lsi",
            verbose = FALSE
        )
        
        # Perform Harmony integration
        co_embed_list[[i]] <- RunHarmony(
            co_embed_list[[i]],
            group.by.vars = "tech",
            reduction.use = "pca",
            max.iter.harmony = 30,
            assay.use = "RNA",
            dims.use = 1:30,
            project.dim = FALSE,
            plot_convergence = FALSE
        )
        
        # Convert assays to "Assay5" after harmonization
        co_embed_list[[i]][["GeneActivity"]] <- as(co_embed_list[[i]][["GeneActivity"]], "Assay5")
        co_embed_list[[i]][["integrated"]] <- as(co_embed_list[[i]][["integrated"]], "Assay5")
        co_embed_list[[i]][["peaks"]] <- as(co_embed_list[[i]][["peaks"]], "Assay5")
    }
    
    # Return the co-embedded list
    return(co_embed_list)
}

#' @export
process_co_embed_list <- function(co_embed_list, RNA_objs, ATAC_objs) {
  pair_list <- list()
  final_pair_list <- list()

  for (i in 1:length(co_embed_list)) {
    if ("ATAC" %in% names(co_embed_list[[i]]@assays)) {
      message(paste0("Processing co-embedded object ", i, ": 'ATAC' assay found."))
    } else if ("peaks" %in% names(co_embed_list[[i]]@assays)) {
      message(paste0("Processing co-embedded object ", i, ": 'ATAC' assay missing, using 'peaks' instead."))
      co_embed_list[[i]][["ATAC"]] <- co_embed_list[[i]][["peaks"]]
    } else {
      message(paste0("Skipping co-embedded object ", i, ": neither 'ATAC' nor 'peaks' assays found."))
      next
    }

    message("Pairing")
    pair_list[[i]] <- PairCells(
      object = co_embed_list[[i]],
      reduction = "harmony",
      pair.by = "tech",
      ident1 = "ATAC",
      ident2 = "RNA"
    )

    message("Selecting cells and filtering")
    rna_cells <- pair_list[[i]]$RNA
    atac_cells <- pair_list[[i]]$ATAC

    obj_rna_subset <- subset(RNA_objs[[i]], cells = rna_cells)
    obj_atac_subset <- subset(ATAC_objs[[i]], cells = atac_cells)

    sel_cells <- unique(c(rna_cells, atac_cells))
    co_embed_list[[i]] <- subset(co_embed_list[[i]], cells = sel_cells)

    message("Creating final paired object")
    final_pair_list[[i]] <- CreatePairedObject(
      pair_list[[i]],
      obj.coembed = co_embed_list[[i]],
      obj.rna = obj_rna_subset,
      obj.atac = obj_atac_subset,
      rna.assay = "RNA",
      atac.assay = "peaks"
    )
  }

  if (length(final_pair_list) == 0) {
    stop("No valid paired objects were created.")
  }

  all_coembed_merge <- merge(final_pair_list[[1]], y = final_pair_list[-1])
  all_coembed_merge <- NormalizeData(all_coembed_merge)
  all_coembed_merge <- all_coembed_merge %>%
    ScaleData() %>%
    FindVariableFeatures() %>%
    RunPCA(verbose = FALSE) %>%
    FindNeighbors(dims = 1:30, verbose = FALSE) %>%
    RunUMAP(dims = 2:25, verbose = FALSE)

  return(all_coembed_merge)
}
