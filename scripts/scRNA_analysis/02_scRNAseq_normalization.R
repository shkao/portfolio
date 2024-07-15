pacman::p_load(scater, scran, sctransform, tidyverse, BiocParallel, patchwork)
setwd("~/Downloads/GSE202109_RAW/")

# Set up parallel processing parameters
bp_params <- MulticoreParam(workers = multicoreWorkers())

# Filter and prepare sample data for plotting
sample_data <- colData(sce) %>%
  as.data.frame() %>%
  filter(Sample == "GSM6094652") %>%
  select(Sample, Barcode, sum) %>%
  mutate(cell_index = 1:n())

# Plot UMI counts before normalization
plot_before_norm <- ggplot(data = sample_data, aes(x = cell_index, y = sum)) +
  geom_bar(stat = 'identity') +
  labs(x = 'Cell Index', y = 'Cell UMI counts', title = "PBMMC_1: Before Normalization") +
  theme_classic() +
  theme(plot.title = element_text(
    hjust = 0.5,
    size = 20,
    color = 'red'
  ))

plot_before_norm

# Normalize UMI counts distribution

# Perform deconvolution to identify clusters
set.seed(100)
clusters <- quickCluster(sce, BPPARAM = bp_params)
table(clusters)

# Compute size factors using deconvolution method
sce <- computePooledFactors(sce,
                            clusters = clusters,
                            min.mean = 0.1,
                            BPPARAM = bp_params)
deconv_size_factors <- sizeFactors(sce)
summary(deconv_size_factors)

# Compute library size factors
lib_size_factors <- librarySizeFactors(sce)
size_factors_df <- data.frame(
  LibrarySizeFactors = lib_size_factors,
  DeconvolutionSizeFactors = deconv_size_factors,
  SampleGroup = sce$Group
)

# Plot comparison of library size factors and deconvolution size factors
ggplot(size_factors_df,
       aes(x = LibrarySizeFactors, y = DeconvolutionSizeFactors)) +
  geom_point(aes(col = SampleGroup)) +
  geom_abline(slope = 1, intercept = 0)

# Apply size factors and normalize counts
sce <- logNormCounts(sce)
assayNames(sce)

# Extract and summarize normalized counts
normalized_counts <- logNormCounts(sce, transform = 'none') %>%
  assay('normcounts') %>%
  as.matrix() %>%
  colSums()

# Prepare data for plotting normalized counts
normalized_counts_df <- tibble(Barcode = names(normalized_counts),
                               normCounts = log2(normalized_counts)) %>%
  inner_join(sample_data, by = 'Barcode')

# Plot UMI counts after normalization
plot_after_norm <- ggplot(data = normalized_counts_df, aes(x = cell_index, y = normCounts)) +
  geom_bar(stat = 'identity') +
  labs(x = 'Cell Index', y = 'Normalized Cell UMI counts', title = "PBMMC_1: After Normalization") +
  theme_classic() +
  theme(plot.title = element_text(
    hjust = 0.5,
    size = 20,
    color = 'red'
  ))

plot_before_norm + plot_after_norm

# Apply sctransform for variance stabilization
umi_counts <- counts(sce)
class(umi_counts)

# Compute gene attributes
gene_attributes <- data.frame(
  mean = rowMeans(umi_counts),
  detection_rate = rowMeans(umi_counts > 0),
  variance = rowVars(umi_counts)
) %>%
  mutate(log_mean = log10(mean),
         log_variance = log10(variance))

dim(gene_attributes)
head(gene_attributes)

# Compute cell attributes
cell_attributes <- data.frame(total_umi = colSums(umi_counts),
                              total_genes = colSums(umi_counts > 0))

dim(cell_attributes)
head(cell_attributes)

# Plot mean-variance relationship for genes
ggplot(gene_attributes, aes(log_mean, log_variance)) +
  geom_point(alpha = 0.3, shape = 16) +
  geom_abline(intercept = 0,
              slope = 1,
              color = 'red')

# Plot mean-detection-rate relationship for genes
x_seq <- seq(from = -3,
             to = 2,
             length.out = 1000)
poisson_model <- data.frame(log_mean = x_seq,
                            detection_rate = 1 - dpois(0, lambda = 10 ^
                                                         x_seq))

ggplot(gene_attributes, aes(log_mean, detection_rate)) +
  geom_point(alpha = 0.3, shape = 16) +
  geom_line(data = poisson_model, color = 'red') +
  theme_gray(base_size = 8)

# Plot cell attributes: total UMI counts vs total genes detected
ggplot(cell_attributes, aes(total_umi, total_genes)) +
  geom_point(alpha = 0.3, shape = 16) +
  geom_density_2d(size = 0.3)

# Estimate model parameters and transform data using sctransform
set.seed(44)
options(future.globals.maxSize = 838860800)

vst_output <- vst(
  umi = umi_counts,
  latent_var = c('log_umi'),
  return_gene_attr = TRUE,
  return_cell_attr = TRUE
)
