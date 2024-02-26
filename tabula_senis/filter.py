import anndata
import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np


def filter_data(adata, gene_threshold=3, cell_threshold=250):
    X_filter = (adata.X > 0).astype(np.int64)
    keep_gene_idx = np.array((X_filter.sum(axis=0) > gene_threshold)).flatten()
    keep_cell_idx = np.array((X_filter.sum(axis=1) > cell_threshold)).flatten()
    return adata[:, keep_gene_idx][keep_cell_idx, :]


def filter_ctype(adata, cell_type):
    task_label_arr = np.array(pd.Series(adata.obs.cell_type.unique())[
                 pd.Series(adata.obs.cell_type.unique()).str.contains(f'{cell_type}')].tolist())
    adata_sub = adata[adata.obs.cell_type.isin(task_label_arr)]
    return adata_sub


def filter_tissue(adata_sub, tissue):
    return adata_sub[adata_sub.obs.tissue.isin([f'{tissue}'])]


def filter_mean_cols(log_mean_df, n_zeros=2):
    zero_mask = (log_mean_df == 0)
    # Count the number of zero entries per column
    zero_counts = zero_mask.sum(axis=0)
    mean_filter_df = log_mean_df.loc[:, zero_counts <= n_zeros]
    return mean_filter_df


def filter_coeff_var(cv_df):
    cv_filter_df = cv_df.loc[:, cv_df.sum(axis=0) > 0]
    df_filtered = cv_filter_df.loc[:, (cv_filter_df != 0).all(axis=0)]
    return df_filtered
