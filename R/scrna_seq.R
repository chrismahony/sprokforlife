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

  anchors <- FindIntegrationAnchors(object.list = scrna.list, dims = 1:50)
  aggr <- IntegrateData(anchorset = anchors, dims = 1:50)

  aggr <- FindVariableFeatures(aggr)
  aggr <- ScaleData(aggr, verbose = FALSE)
  aggr <- RunPCA(aggr, verbose = FALSE)
  aggr <- RunUMAP(aggr, dims = 1:n_dims)
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
  assign("aggr", aggr, envir = .GlobalEnv)
  assign("markers", markers, envir = .GlobalEnv)

  return(list(aggr = aggr, markers = markers))
}




#' @export
process_scATAC_data <- function(h5_file, 
                                metadata_file, 
                                fragments_file, 
                                libraries_map_file, 
                                annotation_db = EnsDb.Hsapiens.v86, 
                                min_cells = 10, 
                                min_features = 200, 
                                tss_threshold = 3, 
                                nucleosome_threshold = 4, 
                                ncount_min = 3000, 
                                ncount_max = 50000, 
                                pct_reads_min = 15, 
                                tss_min = 2, 
                                resolutions = c(0.01, 0.05, 0.1, 0.2)) {
  
  library(Seurat)
  library(Signac)
  library(EnsDb.Hsapiens.v86)
  library(readr)
  
  # Load counts and metadata
  counts <- Read10X_h5(filename = h5_file)
  metadata <- read.csv(file = metadata_file, header = TRUE, row.names = 1)
  
  # Get gene annotations
  annotation <- GetGRangesFromEnsDb(ensdb = annotation_db)
  seqlevels(annotation) <- paste0('chr', seqlevels(annotation))
  
  # Create ChromatinAssay
  chrom_assay <- CreateChromatinAssay(
    counts = counts,
    sep = c(":", "-"),
    fragments = fragments_file,
    min.cells = min_cells,
    min.features = min_features,
    annotation = annotation
  )
  
  # Create Seurat object
  aggr <- CreateSeuratObject(
    counts = chrom_assay,
    assay = "peaks",
    meta.data = metadata
  )
  
  # Add sample names
  cellcodes <- data.frame(barcodes = colnames(aggr))
  cellcodes$libcodes <- gsub(pattern = ".+-", replacement = "", cellcodes$barcodes)
  libraries_map <- read_csv(libraries_map_file)
  cellcodes$samples <- libraries_map$library_id[match(cellcodes$libcodes, libraries_map$library_id)]
  cellcodes <- cellcodes[, "samples", drop = FALSE]
  colnames(cellcodes) <- "orig.ident"
  aggr <- AddMetaData(aggr, cellcodes)
  
  # Compute quality control metrics
  aggr <- NucleosomeSignal(object = aggr)
  aggr <- TSSEnrichment(object = aggr, fast = TRUE)
  aggr$pct_reads_in_peaks <- aggr$peak_region_fragments / aggr$passed_filters * 100
  aggr$high.tss <- ifelse(aggr$TSS.enrichment > tss_threshold, 'High', 'Low')
  aggr$nucleosome_group <- ifelse(aggr$nucleosome_signal > nucleosome_threshold, 'NS > 4', 'NS < 4')
  
  # Filter cells
  aggr <- subset(
    x = aggr,
    subset = nCount_peaks > ncount_min &
      nCount_peaks < ncount_max &
      pct_reads_in_peaks > pct_reads_min &
      nucleosome_signal < nucleosome_threshold &
      TSS.enrichment > tss_min
  )
  
  # Dimensionality reduction and clustering
  aggr <- RunTFIDF(aggr)
  aggr <- FindTopFeatures(aggr, min.cutoff = 'q5')
  aggr <- RunSVD(aggr)
  aggr <- RunUMAP(object = aggr, reduction = 'lsi', dims = 2:30)
  aggr <- FindNeighbors(object = aggr, reduction = 'lsi', dims = 2:30)
  aggr <- FindClusters(object = aggr, verbose = FALSE, algorithm = 3, resolution = resolutions)
  
  # Gene activity analysis
  gene.activities <- GeneActivity(aggr)
  aggr[['RNA']] <- CreateAssayObject(counts = gene.activities)
  aggr <- NormalizeData(
    object = aggr,
    assay = 'RNA',
    normalization.method = 'LogNormalize',
    scale.factor = median(aggr$nCount_RNA)
  )
  
  return(aggr)
}
