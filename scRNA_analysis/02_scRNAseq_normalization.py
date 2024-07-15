import numpy as np
import anndata as ad
import scanpy as sc


def log_normalize(data: ad.AnnData, scale_factor: int) -> ad.AnnData:
    """Normalizes data using a log transformation.

    Args:
        data (ad.AnnData): Dataset
        scale_factor (int): Scale factor for cell-level normalization.

    Returns:
        ad.AnnData: Normalized dataset with normalized counts in data.X
                    and original counts in data.layers["counts"].
    """

    # Recalculate count sum for each cell
    count_sums = np.sum(data.X, axis=1)
    data.obs["n_counts_recalc"] = count_sums

    # Save original counts in a separate layer
    data.layers["counts"] = data.X.copy()

    # Normalize and log transform counts
    data.X = (data.X / count_sums[:, None]) * scale_factor
    sc.pp.log1p(data)

    return data
