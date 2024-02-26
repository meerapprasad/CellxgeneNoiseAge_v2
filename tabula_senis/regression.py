import anndata
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import numpy as np

from tabula_senis.filter import filter_mean_cols
from tabula_senis.plots import plot_linear_regression


def compute_log_mean_and_var(adata_tissue, sorted_ages_mm):
    # fano factor: var/ mean
    log_mean_df = pd.DataFrame(np.zeros((len(adata_tissue.obs.age.unique().tolist()), len(adata_tissue.var_names))),
                               index=sorted_ages_mm, columns=adata_tissue.var['feature_name'])
    for val in sorted_ages_mm:
        log_mean_df.loc[val, :] = np.log1p(adata_tissue[adata_tissue.obs.age == val].X.todense()).mean(axis=0)

    log_var_df = pd.DataFrame(np.zeros((len(adata_tissue.obs.age.unique().tolist()), len(adata_tissue.var_names))),
                              index=sorted_ages_mm, columns=adata_tissue.var['feature_name'])
    for val in sorted_ages_mm:
        log_data = np.log1p(adata_tissue[adata_tissue.obs.age == val].X.todense())
        # Compute squared deviations from the mean
        squared_deviations = (np.array(log_data) - np.array(log_mean_df.loc[val, :])) ** 2
        # Sum squared deviations and divide by the number of cells to get variance
        variance = np.sum(squared_deviations, axis=0) / (log_data.shape[0] - 1)
        # Assign variance to the corresponding row in log_var_df
        log_var_df.loc[val, :] = variance

    return log_mean_df, log_var_df


def linear_regression_on_means(adata_tissue, sorted_ages, output_dir):
    log_mean_df, _ = compute_log_mean_and_var(adata_tissue, sorted_ages)
    # Filter columns with more than 2 zero entries
    mean_filter_df = filter_mean_cols(log_mean_df)
    # binarize the df
    bin_df = (mean_filter_df > 0).astype(int)
    time_points = np.array(mean_filter_df.index.str.extract('(\d+)m')[0], dtype=int)
    time_mean = np.mean(time_points[:, np.newaxis] * bin_df, axis=0)
    # reshape counts_mean to be the same shape as mean_filter_df
    counts_mean = np.mean(mean_filter_df.values, axis=0) * bin_df
    time_val = time_points[:, np.newaxis] * bin_df.values - time_mean.values * bin_df.values
    mean_val = mean_filter_df - counts_mean * bin_df.values
    # Calculate slope (m), for genes without data at an age, it will not affect the mean
    slope = np.sum(time_val * mean_val.values, axis=0) / np.sum((time_val) ** 2, axis=0)

    # Calculate intercept (b)
    intercept = np.mean(mean_filter_df.values, axis=0) - slope * time_mean
    # mean_filter_df.iloc[:, np.argwhere(slope == slope.min())[0][0]]
    # todo: plot the top 100 genes with highest and lowest slopes, and the curve that is fitted
    plot_linear_regression(slope, intercept, mean_filter_df, output_dir, n_plots=100)
    return slope, intercept, mean_filter_df


def linear_regression_least_squares(adata_tissue, sorted_ages, output_dir):
    log_mean_df, _ = compute_log_mean_and_var(adata_tissue, sorted_ages)
    # Filter columns with more than 2 zero entries
    mean_filter_df = filter_mean_cols(log_mean_df)
    # Number of data points
    bin_df = (mean_filter_df > 0).astype(int)
    n = bin_df.sum(axis=0)
    time_points = np.array(mean_filter_df.index.str.extract('(\d+)m')[0], dtype=int)
    # Calculating sums and sum of products
    x = time_points[:, np.newaxis] * bin_df
    y = mean_filter_df.values * bin_df
    sum_x = np.sum(x, axis=0)
    sum_y = np.sum(y, axis=0)
    sum_xy = np.sum(x * y)
    sum_x2 = np.sum(x ** 2)

    # Calculating slope (m) and y-intercept (b)
    slope = (n * sum_xy - sum_x * sum_y) / (n * sum_x2 - sum_x ** 2)
    intercept = (sum_y - slope * sum_x) / n
    return slope, intercept, mean_filter_df