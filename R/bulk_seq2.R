#' @export
run_deseq_contrast_heatmap <- function(dds, meta_data, targetvar, p_value_threshold = 0.05,
                                       log2fc_threshold = 1, condition_colors = NULL,
                                       condition_order = NULL, go_analysis = TRUE,
                                       organism = "mm", annotation_gs) {

  # Generate pairwise comparisons for target variable
  conds <- data.frame(t(combn(unique(as.character(meta_data[[targetvar]])), 2)))

  # Apply DESeq2 on each comparison
  results_list <- lapply(seq_len(nrow(conds)), function(i) {
    cp <- conds[i, ]
    message("Processing comparison: ", cp[1], " vs ", cp[2])
    res <- as.data.frame(DESeq2::results(dds, contrast = c(targetvar, cp[1], cp[2])))
    res[["gene"]] <- rownames(res)
    res[["comparison"]] <- paste0(cp[1], "_vs_", cp[2])
    return(res)
  })
  all_results <- do.call(rbind, results_list)

  # Filter significant results and add a score
  filtered_res <- all_results %>%
    filter(padj < p_value_threshold, abs(log2FoldChange) > log2fc_threshold) %>%
    mutate(score = log2FoldChange * -log10(pvalue)) %>%
    arrange(desc(abs(score)))

  assign("filtered_res", filtered_res, envir = .GlobalEnv)

  # Add clustering information based on log2FoldChange
  filtered_res$cluster <- ifelse(filtered_res$log2FoldChange > 0.001, "UP",
                                 ifelse(filtered_res$log2FoldChange < -0.001, "DOWN", "NO"))

  if (go_analysis) {
    # Prepare gene IDs for GO analysis
    FilteredGeneID <- unique(all_results$gene)
    index <- match(FilteredGeneID, annotation_gs$gene_name)
    ensemblUni <- na.omit(annotation_gs$ensembl_id[index])

    # Match gene names and filter for GO analysis
    index <- match(filtered_res$gene, annotation_gs$gene_name)
    filtered_res$ensembl <- annotation_gs$ensembl_id[index]
    filtered_res <- na.omit(filtered_res)
    filtered_res_pos <- filtered_res %>% filter(log2FoldChange > 0.001)

    # Run GO analysis and plot top results
    go.results <- runGO.all(
      results = filtered_res_pos,
      background_ids = ensemblUni,
      gene_id_col = "ensembl",
      gene_id_type = "ensembl",
      sample_col = "comparison",
      p_col = "padj",
      p_threshold = 0.05,
      species = organism
    )
    go.results <- filterGenesets(go.results)
    go.results.top <- go.results %>%
      group_by(comparison) %>%
      top_n(n = 5, -p.val)

    sampleEnrichmentDotplot(go.results.top,
                            selection_col = "description",
                            selected_genesets = unique(go.results.top$description),
                            sample_id_col = "comparison",
                            fill_var = "odds.ratio",
                            maxl = 50,
                            title = "GO Term Enrichment",
                            rotate_sample_labels = TRUE)

    assign("go.results", go.results, envir = .GlobalEnv)
    assign("go.results.top", go.results.top, envir = .GlobalEnv)
  }

  # Generate VST-transformed data and create heatmap
  vsd <- vst(dds, blind = TRUE)
  gene_filter <- unique(filtered_res$gene)
  vsd_mat <- assay(vsd)[rownames(vsd) %in% gene_filter, ]
  scale_vsd <- t(scale(t(vsd_mat)))

  # Reorder conditions if specified
  if (!is.null(condition_order)) {
    meta_data[[targetvar]] <- factor(meta_data[[targetvar]], levels = condition_order)
    ordered_indices <- order(meta_data[[targetvar]])
    scale_vsd <- scale_vsd[, ordered_indices]
    meta_data <- meta_data[ordered_indices, ]
  }

  # Default colors for condition annotations if not provided
  if (is.null(condition_colors)) {
    condition_colors <- RColorBrewer::brewer.pal(n = length(unique(meta_data[[targetvar]])), name = "Set1")
    names(condition_colors) <- unique(meta_data[[targetvar]])
  }

  top_annotation <- HeatmapAnnotation(
    condition = meta_data[[targetvar]],
    col = list(condition = condition_colors)
  )

  Heatmap(scale_vsd,
          top_annotation = top_annotation,
          col = colorRamp2(c(-1.5, 0, 1.5), c("blue", "white", "red")),
          row_names_gp = gpar(fontsize = 4),
          cluster_columns = FALSE,
          cluster_rows = TRUE,
          show_row_names = FALSE,
          show_column_names = FALSE,
          border = TRUE)
}

