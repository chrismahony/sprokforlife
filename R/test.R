#' @export
process_deseq_PCA_NEW <- function(counts, meta_data, relevel_condition = NULL, min_gene_count = 200) {
  # Check if batch column is present in meta_data
  if ("batch" %in% colnames(meta_data)) {
    message("Batch correction using ComBat_seq...")
    library(sva)
    batch <- meta_data$batch
    adjusted <- ComBat_seq(counts, batch = batch, group = NULL)
  } else {
    message("No batch column detected, skipping batch correction...")
    adjusted <- counts
  }

  # Create DESeq2 dataset
  dds <- DESeqDataSetFromMatrix(countData = adjusted,
                                colData = meta_data,
                                design = ~ 1)

  # Filter lowly expressed genes based on min_gene_count
  message(paste("Filtering genes with counts below", min_gene_count, "..."))
  keep <- rowSums(counts(dds)) > min_gene_count
  dds <- dds[keep, ]

  # Relevel condition column if a new order is provided
  if (!is.null(relevel_condition)) {
    message("Releveling condition factor to specified order...")
    if (!all(relevel_condition %in% levels(meta_data$condition))) {
      stop("The specified relevel order contains values not present in the 'condition' column.")
    }
    dds@colData[['condition']] <- factor(dds@colData[['condition']], levels = relevel_condition)
  }

  # Set design formula and run DESeq
  design(dds) <- formula(~ condition)

  message("Running DESeq analysis...")
  dds <- DESeq(dds, test = "Wald")

  # Plot dispersion estimates
  plotDispEsts(dds)

  # PCA analysis
  message("Performing PCA analysis...")
  rld_v3 <- rlog(dds, blind = TRUE)

  # Plot PCA using ggplot2
  pcaData <- plotPCA(rld_v3, intgroup = c("condition"), returnData = TRUE)
  percentVar <- round(100 * attr(pcaData, "percentVar"))

  p <- ggplot(pcaData, aes(PC1, PC2, color = condition)) +
    geom_point(size = 3) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    coord_fixed()

  print(p)

  return(list(dds = dds, pca_plot = p))
}
