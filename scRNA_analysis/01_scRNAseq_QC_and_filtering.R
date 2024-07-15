pacman::p_load(DropletUtils,
               scater,
               ensembldb,
               AnnotationHub,
               BiocParallel,
               tidyverse,
               stringr)
setwd("~/Downloads/GSE202109_RAW/")

# Load sample metadata
sample_metadata <- read.csv("GSE202109_metadata.csv.gz")

bp_params <- MulticoreParam(workers = multicoreWorkers())
file_list <- dir(pattern = "GSM")
names(file_list) <- str_split_fixed(file_list, "_", 2)[, 1]
sce <- read10xCounts(file_list, col.names = TRUE, BPPARAM = bp_params)

# Check samples in the dataset
colData(sce) %>%
  as.data.frame() %>%
  select(Sample) %>%
  distinct()

sce

# Modify the droplet annotation
sce$Barcode <- rownames(colData(sce))
colData(sce) <- merge(colData(sce), sample_metadata, by = "Sample", sort = FALSE)
rownames(colData(sce)) <- sce$Barcode

# Remove undetected genes
detected_genes <- rowSums(counts(sce)) > 0
sce <- sce[detected_genes, ]
mean(detected_genes)
# 0.835 -> approximately 83.5% of genes were detected in at least one sample

# Annotate genes
annotation_hub <- AnnotationHub()
query_results <- query(annotation_hub, c("Homo sapiens", "EnsDb", 108))[[1]]

gene_ids <- rowData(sce)$ID
gene_annotations <- AnnotationDbi::select(
  query_results,
  keys = gene_ids,
  keytype = "GENEID",
  columns = c("GENEID", "SEQNAME")
) %>%
  set_names(c("ID", "Chromosome"))

rowData(sce) <- merge(
  rowData(sce),
  gene_annotations,
  by = "ID",
  sort = FALSE,
  all.x = TRUE
)
rownames(rowData(sce)) <- rowData(sce)$ID

rowData(sce)

# Add per cell QC metrics
mitochondrial_genes <- which(rowData(sce)$Chromosome == "MT")
sce <- addPerCellQC(sce, subsets = list(Mito = mitochondrial_genes))

plotColData(sce,
            x = "Sample",
            y = "sum",
            other_fields = "Group") +
  facet_wrap( ~ Group, nrow = 1, scales = "free_x") +
  scale_y_log10() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  ggtitle("Total count")

plotColData(sce,
            x = "Sample",
            y = "detected",
            other_fields = "Group") +
  facet_wrap( ~ Group, nrow = 1, scales = "free_x") +
  scale_y_log10() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  ggtitle("Detected features")

plotColData(sce,
            x = "Sample",
            y = "subsets_Mito_percent",
            other_fields = "Group") +
  facet_wrap( ~ Group, nrow = 1, scales = "free_x") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  ggtitle("Mito percent")

colData(sce) %>%
  as.data.frame() %>%
  arrange(subsets_Mito_percent) %>%
  ggplot(aes(x = sum, y = detected)) +
  geom_point(aes(colour = subsets_Mito_percent > 10)) +
  facet_wrap(vars(Group))

# Filter low-quality cells based on three metrics: library size, number of features, and mitochondrial percent
cell_qc_filters <- quickPerCellQC(
  colData(sce),
  percent_subsets = c("subsets_Mito_percent"),
  batch = sce$Sample
)

as.data.frame(cell_qc_filters) %>%
  summarise(across(everything(), sum))

# Add the columns in the droplet annotation with these new filters
colData(sce) <- cbind(colData(sce), cell_qc_filters)

plotColData(
  sce,
  x = "Sample",
  y = "sum",
  other_fields = "Group",
  colour_by = "low_lib_size"
) +
  facet_wrap(vars(Group), nrow = 1, scales = "free_x") +
  scale_y_log10() +
  labs(y = "Total count", title = "Total count") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  guides(colour = guide_legend(title = "Discarded"))

plotColData(
  sce,
  x = "Sample",
  y = "detected",
  other_fields = "Group",
  colour_by = "low_n_features"
) +
  facet_wrap(vars(Group), nrow = 1, scales = "free_x") +
  scale_y_log10() +
  labs(y = "Genes detected", title = "Genes detected") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  guides(colour = guide_legend(title = "Discarded"))

plotColData(
  sce,
  x = "Sample",
  y = "subsets_Mito_percent",
  other_fields = "Group",
  colour_by = "high_subsets_Mito_percent"
) +
  facet_wrap(vars(Group), nrow = 1, scales = "free_x") +
  labs(y = "Percentage mitochondrial UMIs", title = "Mitochondrial UMIs") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  guides(colour = guide_legend(title = "Discarded"))

# Filter out poor quality cells
sce <- sce[, !sce$discard]
