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


def filter_tissue(adata, tissue):
    return adata[adata.obs.tissue.isin([f'{tissue}'])]

# todo: ensure that first entry is not zero
def filter_mean_cols(log_mean_df, n_zeros=1):
    zero_mask = (log_mean_df == 0)
    # Count the number of zero entries per column
    zero_counts = zero_mask.sum(axis=0)
    mean_filter_df = log_mean_df.loc[:, zero_counts <= n_zeros]
    mean_filter_df2 = mean_filter_df.loc[:,~(mean_filter_df.iloc[0] == 0)]
    return mean_filter_df2


def filter_coeff_var(cv_df):
    cv_filter_df = cv_df.loc[:, cv_df.sum(axis=0) > 0]
    df_filtered = cv_filter_df.loc[:, (cv_filter_df != 0).all(axis=0)]
    return df_filtered


def filter_age(adata_sub):
    return adata_sub[adata_sub.obs.age.isin(['3m', '18m', '21m', '24m', '30m'])]


def combine_labels_adata(adata_sub, age_map_dict, age_col):
    adata_sub.obs[age_col] = adata_sub.obs['age'].map(age_map_dict)
    return adata_sub


def extract_names(mean_filter_df):
    extracted = mean_filter_df.index.str.extract(r'(\d+)-?(\d+)?')
    # Convert extracted strings to numeric (float for NaN compatibility)
    extracted = extracted.apply(pd.to_numeric)
    # If there's a second number (i.e., a range), calculate the average, otherwise use the first number
    time_points = np.where(extracted[1].isna(), extracted[0], (extracted[0] + extracted[1]) / 2)
    return time_points

