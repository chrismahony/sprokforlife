#' @export
run_deseq_contrast_heatmap <- function(dds, meta_data, targetvar, p_value_threshold = 0.05,
                                       log2fc_threshold = 1, condition_colors = NULL,
                                       condition_order = NULL, go_analysis = TRUE,
                                       organism = "mm", annotation_gs) {

  # Generate all pairwise comparisons for the target variable
  comps1 <- data.frame(t(combn(unique(as.character(meta_data[[targetvar]])), 2)))

  # Perform DESeq2 analysis for each comparison and collect results
  ress <- apply(comps1, 1, function(cp) {
    print(cp)
    res <- data.frame(DESeq2::results(dds, contrast = c(targetvar, cp[1], cp[2])))
    res[["gene"]] <- rownames(res)
    res[["comparison"]] <- paste0(cp[1], "_vs_", cp[2])
    res
  })

  # Combine all results into a single dataframe
  res1 <- Reduce(rbind, ress)

  # Filter significant results and calculate a "score" based on log2FC and p-value
  subres <- res1 %>%
    filter(padj < p_value_threshold, abs(log2FoldChange) > log2fc_threshold) %>%
    mutate(score = log2FoldChange * (-log10(pvalue))) %>%
    arrange(desc(abs(score)))

  assign("subres", subres, envir = .GlobalEnv)

  # Add cluster assignments based on log2FoldChange
  subres$cluster <- "NO"
  subres$cluster[subres$log2FoldChange > 0.001] <- "UP"
  subres$cluster[subres$log2FoldChange < -0.001] <- "DOWN"

  # Perform GO analysis if specified
  go.results <- NULL  # Initialize variable
  go.results.top <- NULL  # Initialize variable

  assign("res1", res1, envir = .GlobalEnv)


  if (go_analysis) {
    # Prepare the unique gene IDs for background
    FilteredGeneID <- unique(res1$gene)
    index <- match(FilteredGeneID, annotation_gs$gene_name)
    ensemblUni <- annotation_gs$ensembl_id[index]
    ensemblUni <- na.omit(ensemblUni)

    assign("subres", subres, envir = .GlobalEnv)

    # Match the gene names to obtain Ensembl IDs
    index <- match(subres$gene, annotation_gs$gene_name)
    subres$ensembl <- annotation_gs$ensembl_id[index]
    subres <- na.omit(subres)

    subresf <- subres %>% filter(log2FoldChange >0.001)

    # Run GO analysis
    go.results <- runGO.all(
      results = subresf,
      background_ids = ensemblUni,
      gene_id_col = "ensembl",
      gene_id_type = "ensembl",
      sample_col = "comparison",
      p_col = "padj",
      p_threshold = 0.05,
      species = "mm"
    )

    # Filter genesets and get top results
    go.results <- filterGenesets(go.results)
    go.results.top <- go.results %>% group_by("comparison") %>% top_n(n = 5, -p.val)

    # Create a dot plot for GO term enrichment
    sampleEnrichmentDotplot(go.results.top,
                            selection_col = "description",
                            selected_genesets = unique(go.results.top$description),
                            sample_id_col = "comparison",
                            fill_var = "odds.ratio",
                            maxl = 50,
                            title = "GO Term Enrichment",
                            rotate_sample_labels = TRUE)
  }

  # Save the results into the R environment

  assign("go.results", go.results, envir = .GlobalEnv)
  assign("go.results.top", go.results.top, envir = .GlobalEnv)

  # Generate a heatmap if more than 10 unique genes are found
  if (length(unique(subres$gene)) > 10) {
    # Perform variance stabilizing transformation (VST)
    vsd <- tryCatch({
      vst(dds, blind = TRUE)
    }, error = function(e) {
      message(e)
      print(e)
      return(NULL)
    })

    if (!is.null(vsd)) {
      # Extract the transformed counts matrix
      vsd_mat <- assay(vsd)

      # Filter for the relevant genes (features) from the significant results
      feats <- unique(subres$gene)

      # Subset the matrix to only include the significant genes
      sub_vsd_mat <- vsd_mat[rownames(vsd_mat) %in% feats, ]

      # Scale the subset matrix by gene (row) to have a mean of 0 and variance of 1
      scale_sub_vsd <- t(scale(t(sub_vsd_mat)))

      # Set default colors for condition annotations if not provided
      if (is.null(condition_colors)) {
        condition_colors <- RColorBrewer::brewer.pal(n = length(unique(meta_data[[targetvar]])), name = "Set1")
        names(condition_colors) <- unique(meta_data[[targetvar]])
      }

      # Reorder conditions based on user input
      if (!is.null(condition_order)) {
        # Ensure the levels are reordered based on the input order
        meta_data[[targetvar]] <- factor(meta_data[[targetvar]], levels = condition_order)

        # Get the order of the columns based on the new factor levels
        ordered_indices <- order(meta_data[[targetvar]])

        # Reorder the scaled matrix based on the ordered indices
        scale_sub_vsd <- scale_sub_vsd[, ordered_indices]

        # Reorder the metadata based on the ordered indices for annotations
        meta_data <- meta_data[ordered_indices, ]
      }

      # Define annotation for the condition
      condition_ann <- meta_data[[targetvar]]

      # Combine both condition and sample annotations in the top annotation
      top_annotation <- HeatmapAnnotation(
        condition = condition_ann,
        col = list(condition = condition_colors)
      )

      assign("scale_sub_vsd", scale_sub_vsd, envir = .GlobalEnv)
      assign("top_annotation", top_annotation, envir = .GlobalEnv)



      # Plot the heatmap with top (condition and sample) annotations
      Heatmap(scale_sub_vsd,
              top_annotation = top_annotation,
              col = colorRamp2(c(-1.5, 0, 1.5), c("blue", "white", "red")),
              row_names_gp = gpar(fontsize = 4),
              cluster_columns = FALSE,
              cluster_rows = TRUE,
              show_row_names = FALSE,
              show_column_names = FALSE,
              border = TRUE)
    }
  } else {
    print("Not enough significant genes for heatmap.")
  }
}
