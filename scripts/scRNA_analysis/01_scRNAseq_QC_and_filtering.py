import numpy as np
import scanpy as sc
from typing import Optional
import anndata as ad


def compute_mitochondrial_fraction(adata: ad.AnnData) -> None:
    """Compute the fraction of mitochondrial genes for each cell and store in adata.obs['mito_pct']."""
    mito_genes = adata.var_names.str.startswith("MT-")
    mito_counts = np.sum(adata[:, mito_genes].X, axis=1).A1
    total_counts = np.sum(adata.X, axis=1).A1
    adata.obs["mito_pct"] = (mito_counts / total_counts) * 100.0


def filter_by_mito_fraction(adata: ad.AnnData, max_mito_pct: float) -> ad.AnnData:
    """Filter cells based on the maximum allowed mitochondrial gene fraction."""
    if "mito_pct" not in adata.obs:
        compute_mitochondrial_fraction(adata)
    return adata[adata.obs["mito_pct"] < max_mito_pct]


def filter_by_counts(
    adata: ad.AnnData, min_counts: Optional[int], max_counts: Optional[int]
) -> ad.AnnData:
    """Filter cells based on the number of counts."""
    adata.obs["total_counts"] = np.sum(adata.X, axis=1).A1
    sc.pp.filter_cells(adata, min_counts=min_counts, max_counts=max_counts)
    return adata


def compute_cells_per_gene(adata: ad.AnnData) -> None:
    """Compute the number of cells expressing each gene and store in adata.var['n_cells']."""
    adata.var["n_cells"] = np.sum(adata.X > 0, axis=0).A1


def filter_genes_and_cells(
    adata: ad.AnnData,
    min_genes: Optional[int],
    max_genes: Optional[int],
    min_cells: Optional[int],
) -> ad.AnnData:
    """Filter cells and genes based on the provided thresholds."""
    adata.obs["n_genes"] = np.sum(adata.X > 0, axis=1).A1
    sc.pp.filter_cells(adata, min_genes=min_genes, max_genes=max_genes)
    if min_cells is not None:
        sc.pp.filter_genes(adata, min_cells=min_cells)
    return adata


def perform_standard_qc(adata: ad.AnnData, qc_params: dict) -> ad.AnnData:
    """Perform standard quality control on the dataset based on the provided parameters."""
    if "max_mito_pct" in qc_params:
        adata = filter_by_mito_fraction(adata, qc_params["max_mito_pct"])
    adata = filter_by_counts(
        adata, qc_params.get("min_counts"), qc_params.get("max_counts")
    )
    adata = filter_genes_and_cells(
        adata,
        qc_params.get("min_genes"),
        qc_params.get("max_genes"),
        qc_params.get("min_cells"),
    )
    return adata
