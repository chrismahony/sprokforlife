#' @export
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
