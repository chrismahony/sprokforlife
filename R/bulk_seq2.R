#' @export
run_deseq_contrast_heatmap  <- function (dds, meta_data, targetvar, p_value_threshold = 0.05, 
    log2fc_threshold = 1, condition_colors = NULL, condition_order = NULL, 
    go_analysis = TRUE, organism = "mm", annotation_gs) 
{
    conds <- data.frame(t(combn(unique(as.character(meta_data[[targetvar]])), 
        2)))
    results_df <- apply(conds, 1, function(cp) {
        print(cp)
        results_df <- data.frame(DESeq2::results(dds, contrast = c(targetvar, 
            cp[2], cp[1])))
        results_df[["gene"]] <- rownames(results_df)
        results_df[["comparison"]] <- paste0(cp[2], "_vs_", cp[1])
        results_df
    })
    results_all <- Reduce(rbind, results_df)
    DEGs_f <- results_all %>% filter(padj < p_value_threshold, 
        abs(log2FoldChange) > log2fc_threshold) %>% mutate(score = log2FoldChange * 
        (-log10(pvalue))) %>% arrange(desc(abs(score)))
    assign("DEGs_f", DEGs_f, envir = .GlobalEnv)
    DEGs_f$cluster <- "NO"
    DEGs_f$cluster[DEGs_f$log2FoldChange > 0.001] <- "UP"
    DEGs_f$cluster[DEGs_f$log2FoldChange < -0.001] <- "DOWN"
    go.results <- NULL
    go.results.top <- NULL
    assign("results_all", results_all, envir = .GlobalEnv)
    if (go_analysis) {
        FilteredGeneID <- unique(results_all$gene)
        index <- match(FilteredGeneID, annotation_gs$gene_name)
        ensemblUni <- annotation_gs$ensembl_id[index]
        ensemblUni <- na.omit(ensemblUni)
        assign("DEGs_f", DEGs_f, envir = .GlobalEnv)
        index <- match(DEGs_f$gene, annotation_gs$gene_name)
        DEGs_f$ensembl <- annotation_gs$ensembl_id[index]
        DEGs_f <- na.omit(DEGs_f)
        DEGs_f_2 <- DEGs_f %>% filter(log2FoldChange > 0.001)
        go.results <- runGO.all(results = DEGs_f_2, background_ids = ensemblUni, 
            gene_id_col = "ensembl", gene_id_type = "ensembl", 
            sample_col = "comparison", p_col = "padj", p_threshold = 0.05, 
            species = organism)
        go.results <- filterGenesets(go.results)
        go.results.top <- go.results %>% group_by(comparison) %>% 
            top_n(n = 5, -p.val)
        sampleEnrichmentDotplot(go.results.top, selection_col = "description", 
            selected_genesets = unique(go.results.top$description), 
            sample_id_col = "comparison", fill_var = "odds.ratio", 
            maxl = 50, title = "GO Term Enrichment", rotate_sample_labels = TRUE)
    }
    assign("go.results", go.results, envir = .GlobalEnv)
    assign("go.results.top", go.results.top, envir = .GlobalEnv)
    vsd <- vst(dds, blind = TRUE)
    gene_filter <- unique(DEGs_f$gene)
    vsd_mat <- assay(vsd)[rownames(vsd) %in% gene_filter, ]
    scale_sub_vsd <- t(scale(t(vsd_mat)))
    if (is.null(condition_colors)) {
        condition_colors <- RColorBrewer::brewer.pal(n = length(unique(meta_data[[targetvar]])), 
            name = "Set1")
        names(condition_colors) <- unique(meta_data[[targetvar]])
    }
    if (!is.null(condition_order)) {
        meta_data[[targetvar]] <- factor(meta_data[[targetvar]], 
            levels = condition_order)
        ordered_indices <- order(meta_data[[targetvar]])
        scale_sub_vsd <- scale_sub_vsd[, ordered_indices]
        meta_data <- meta_data[ordered_indices, ]
    }
    condition_ann <- meta_data[[targetvar]]
    top_annotation <- HeatmapAnnotation(condition = condition_ann, 
        col = list(condition = condition_colors))
    assign("scale_sub_vsd", scale_sub_vsd, envir = .GlobalEnv)
    assign("top_annotation", top_annotation, envir = .GlobalEnv)
   
  print(Heatmap(scale_sub_vsd, top_annotation = top_annotation, col = colorRamp2(c(-1.5, 
        0, 1.5), c("blue", "white", "red")), row_names_gp = gpar(fontsize = 4), 
        cluster_columns = FALSE, cluster_rows = TRUE, show_row_names = FALSE, 
        show_column_names = FALSE, border = TRUE))
  
    sampleEnrichmentDotplot(go.results.top, selection_col = "description", 
        selected_genesets = unique(go.results.top$description), 
        sample_id_col = "comparison", fill_var = "odds.ratio", 
        maxl = 50, title = "Go term", rotate_sample_labels = T)
}




process_scrna_data <- function(file_list, sample_names, n_dims = 50, target_n_clusters = 10,
                               resolution_range = seq(0.05, 0.3, by = 0.5),
                               min_nFeature_RNA = 500, max_nFeature_RNA = 7000,
                               max_percent_mt = 10) {
  data.10x <- list()
  scrna.list <- list()

  if (length(sample_names) != length(file_list)) {
    stop("The length of sample_names must match the length of file_list")
  }

  for (i in 1:length(file_list)) {
    data.10x[[i]] <- Read10X(data.dir = file_list[[i]])
  }

  for (i in 1:length(data.10x)) {
    scrna.list[[i]] <- CreateSeuratObject(counts = data.10x[[i]], min.cells = 3, min.features = 0, project = sample_names[i])
    scrna.list[[i]][["percent.mt"]] <- PercentageFeatureSet(object = scrna.list[[i]], pattern = "^MT-")

    # Subset based on the provided thresholds for nFeature_RNA and percent.mt
    scrna.list[[i]] <- subset(scrna.list[[i]], subset = nFeature_RNA > min_nFeature_RNA &
                                                  nFeature_RNA < max_nFeature_RNA &
                                                  percent.mt < max_percent_mt)

    scrna.list[[i]] <- NormalizeData(object = scrna.list[[i]])
    scrna.list[[i]] <- ScaleData(object = scrna.list[[i]])
    scrna.list[[i]] <- FindVariableFeatures(object = scrna.list[[i]])
    scrna.list[[i]] <- RunPCA(object = scrna.list[[i]], verbose = FALSE)
  }

  names(scrna.list) <- sample_names

  anchors <- FindIntegrationAnchors(object.list = scrna.list, reduction = "rpca", dims = 1:50)
  aggr <- IntegrateData(anchorset = anchors, dims = 1:50)

  aggr <- FindVariableFeatures(aggr)
  aggr <- ScaleData(aggr, verbose = FALSE)
  aggr <- RunPCA(aggr, verbose = FALSE)
  aggr <- RunUMAP(aggr, dims = 1:n_n_dims)
  aggr <- FindNeighbors(aggr, dims = 1:n_dims)

  for (resolution in resolution_range) {
    aggr <- FindClusters(aggr, resolution = resolution, verbose = FALSE)
    num_clusters <- length(unique(Idents(aggr)))

    if (num_clusters >= target_n_clusters) {
      print(paste("Target number of clusters reached at resolution", resolution, "with", num_clusters, "clusters"))
      break
    }
  }


  markers <- FindAllMarkers(aggr, only.pos = TRUE)

  return(list(aggr = aggr, markers = markers))
}
